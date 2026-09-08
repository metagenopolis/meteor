use std::collections::{BTreeMap, HashMap};

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{IndexedReader, Read, Reader, Record};
use rust_htslib::faidx;
use std::fs::File;
use std::io::{BufRead, BufReader};

mod freebayes;
mod vcf;

fn open_cram(cram_path: &str, ref_path: &str) -> PyResult<Reader> {
    let mut reader = Reader::from_path(cram_path)
        .map_err(|e| PyIOError::new_err(format!("failed to open CRAM {cram_path}: {e}")))?;
    reader
        .set_reference(ref_path)
        .map_err(|e| PyIOError::new_err(format!("failed to set CRAM reference {ref_path}: {e}")))?;
    Ok(reader)
}

#[pyfunction]
fn count_cram_records(cram_path: &str, ref_path: &str) -> PyResult<usize> {
    let mut reader = open_cram(cram_path, ref_path)?;
    let mut record = Record::new();
    let mut count: usize = 0;
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| PyIOError::new_err(format!("CRAM read error: {e}")))?;
        count += 1;
    }
    Ok(count)
}

#[pyfunction]
fn sum_cram_cigar_lengths(cram_path: &str, ref_path: &str) -> PyResult<(u64, u64)> {
    let mut reader = open_cram(cram_path, ref_path)?;
    let mut record = Record::new();
    let mut total_len: u64 = 0;
    let mut total_ops: u64 = 0;
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| PyIOError::new_err(format!("CRAM read error: {e}")))?;
        for cigar in record.cigar().iter() {
            total_len += u64::from(cigar.len());
            total_ops += 1;
        }
    }
    Ok((total_len, total_ops))
}

fn cigar_op_to_int(op: Cigar) -> u32 {
    match op {
        Cigar::Match(_) => 0,
        Cigar::Ins(_) => 1,
        Cigar::Del(_) => 2,
        Cigar::RefSkip(_) => 3,
        Cigar::SoftClip(_) => 4,
        Cigar::HardClip(_) => 5,
        Cigar::Pad(_) => 6,
        Cigar::Equal(_) => 7,
        Cigar::Diff(_) => 8,
    }
}

fn cigar_tuples(record: &Record) -> Vec<(u32, u32)> {
    record
        .cigar()
        .iter()
        .map(|cigar| (cigar_op_to_int(*cigar), cigar.len()))
        .collect()
}

fn extract_nm(record: &Record) -> Option<u32> {
    match record.aux(b"NM").ok()? {
        Aux::I8(v) if v >= 0 => Some(v as u32),
        Aux::U8(v) => Some(v as u32),
        Aux::I16(v) if v >= 0 => Some(v as u32),
        Aux::U16(v) => Some(v as u32),
        Aux::I32(v) if v >= 0 => Some(v as u32),
        Aux::U32(v) => Some(v),
        _ => None,
    }
}

fn bytes_to_string(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes)
        .trim_end_matches('\0')
        .to_string()
}

#[pyclass]
#[derive(Clone)]
struct CramRecord {
    #[pyo3(get)]
    query_name: String,
    #[pyo3(get)]
    reference_name: String,
    #[pyo3(get)]
    cigar_tuples: Vec<(u32, u32)>,
    #[pyo3(get)]
    nm: Option<u32>,
}

#[pyfunction]
fn cram_records(cram_path: &str, ref_path: &str) -> PyResult<Vec<CramRecord>> {
    let mut reader = open_cram(cram_path, ref_path)?;
    let header = reader.header().clone();
    let target_names: Vec<String> = header
        .target_names()
        .iter()
        .map(|name| bytes_to_string(name))
        .collect();

    let mut record = Record::new();
    let mut records: Vec<CramRecord> = Vec::new();
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| PyIOError::new_err(format!("CRAM read error: {e}")))?;

        let query_name = bytes_to_string(record.qname());
        let reference_name = if record.tid() >= 0 && (record.tid() as usize) < target_names.len() {
            target_names[record.tid() as usize].clone()
        } else {
            "*".to_string()
        };

        records.push(CramRecord {
            query_name,
            reference_name,
            cigar_tuples: cigar_tuples(&record),
            nm: extract_nm(&record),
        });
    }
    Ok(records)
}

#[pyclass]
#[derive(Clone)]
struct GeneCount {
    #[pyo3(get)]
    gene_id: i32,
    #[pyo3(get)]
    gene_length: i32,
    #[pyo3(get)]
    count: f64,
}

#[pyclass]
#[derive(Clone)]
struct MspCountResult {
    #[pyo3(get)]
    gene_counts: Vec<GeneCount>,
    #[pyo3(get)]
    counted_reads: usize,
}

/// Sum of CIGAR lengths for operators considered "aligned nucleotides" by meteor:
/// M (0), I (1), D (2). RefSkip N (3) is intentionally excluded.
fn aligned_nucleotides(record: &Record) -> u32 {
    record
        .cigar()
        .iter()
        .filter_map(|cigar| match cigar {
            Cigar::Match(len) | Cigar::Ins(len) | Cigar::Del(len) => Some(*len),
            _ => None,
        })
        .sum()
}

#[pyfunction]
fn count_msp(
    cram_path: &str,
    ref_path: &str,
    identity_threshold: f64,
    counting_type: &str,
) -> PyResult<MspCountResult> {
    if !matches!(counting_type, "smart_shared" | "unique" | "total") {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "{counting_type} is not a valid counting type"
        )));
    }

    let mut reader = open_cram(cram_path, ref_path)?;
    let header = reader.header().clone();
    let mut database: BTreeMap<i32, i32> = BTreeMap::new();
    for tid in 0..header.target_count() {
        let name = bytes_to_string(header.target_names()[tid as usize]);
        if let Ok(gene_id) = name.parse::<i32>() {
            let length = header.target_len(tid).unwrap_or(0) as i32;
            database.insert(gene_id, length);
        }
    }

    let mut record = Record::new();
    let mut reads: HashMap<String, (f64, Vec<i32>)> = HashMap::new();

    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| PyIOError::new_err(format!("CRAM read error: {e}")))?;

        let query_name = bytes_to_string(record.qname());
        let reference_name = if record.tid() >= 0 && (record.tid() as u32) < header.target_count() {
            bytes_to_string(header.target_names()[record.tid() as usize])
        } else {
            continue;
        };

        let gene_id = match reference_name.parse::<i32>() {
            Ok(id) => id,
            Err(_) => continue,
        };

        let nm = extract_nm(&record).unwrap_or(0) as f64;
        let aligned = aligned_nucleotides(&record) as f64;
        if aligned <= 0.0 {
            continue;
        }
        let identity = (aligned - nm) / aligned;
        if identity < identity_threshold {
            continue;
        }

        let score = identity;
        match reads.get_mut(&query_name) {
            Some((prev_score, genes)) => {
                if (score - *prev_score).abs() < f64::EPSILON {
                    genes.push(gene_id);
                } else if score > *prev_score {
                    *prev_score = score;
                    *genes = vec![gene_id];
                }
            }
            None => {
                reads.insert(query_name, (score, vec![gene_id]));
            }
        }
    }

    let counted_reads = reads.len();

    if counting_type == "total" {
        let mut abundance: BTreeMap<i32, f64> = database.keys().map(|&g| (g, 0.0)).collect();
        for (_, (_, genes)) in reads {
            for gene in genes {
                *abundance.entry(gene).or_insert(0.0) += 1.0;
            }
        }
        let gene_counts: Vec<GeneCount> = abundance
            .into_iter()
            .map(|(gene_id, count)| GeneCount {
                gene_id,
                gene_length: *database.get(&gene_id).unwrap_or(&0),
                count,
            })
            .collect();
        return Ok(MspCountResult {
            gene_counts,
            counted_reads,
        });
    }

    let mut unique_on_gene: BTreeMap<i32, f64> = database.keys().map(|&g| (g, 0.0)).collect();
    let mut multiple_reads: Vec<(String, Vec<i32>)> = Vec::new();

    for (read_id, (_, genes)) in reads {
        if genes.len() == 1 {
            *unique_on_gene.entry(genes[0]).or_insert(0.0) += 1.0;
        } else {
            multiple_reads.push((read_id, genes));
        }
    }

    if counting_type == "unique" {
        let gene_counts: Vec<GeneCount> = unique_on_gene
            .into_iter()
            .map(|(gene_id, count)| GeneCount {
                gene_id,
                gene_length: *database.get(&gene_id).unwrap_or(&0),
                count,
            })
            .collect();
        return Ok(MspCountResult {
            gene_counts,
            counted_reads,
        });
    }

    let mut co_dict: HashMap<(String, i32), f64> = HashMap::new();
    let mut read_dict: HashMap<i32, Vec<String>> = HashMap::new();

    for (read_id, genes) in &multiple_reads {
        let som: f64 = genes
            .iter()
            .map(|gene| unique_on_gene.get(gene).copied().unwrap_or(0.0))
            .sum();

        if som == 0.0 {
            let n = genes.len() as f64;
            for gene in genes {
                co_dict.insert((read_id.clone(), *gene), 1.0 / n);
                read_dict.entry(*gene).or_default().push(read_id.clone());
            }
            continue;
        }

        let duplicated_genes: std::collections::HashSet<i32> = genes
            .iter()
            .filter(|gene| genes.iter().filter(|g| *g == *gene).count() > 1)
            .copied()
            .collect();

        for gene in genes {
            let nb_unique = unique_on_gene.get(gene).copied().unwrap_or(0.0);
            if nb_unique == 0.0 {
                continue;
            }
            let key = (read_id.clone(), *gene);
            let value = nb_unique / som;
            if duplicated_genes.contains(gene) {
                *co_dict.entry(key).or_insert(0.0) += value;
            } else {
                co_dict.insert(key, value);
            }
            read_dict.entry(*gene).or_default().push(read_id.clone());
        }
    }

    for read_list in read_dict.values_mut() {
        read_list.sort_unstable();
        read_list.dedup();
    }

    let mut abundance = unique_on_gene.clone();
    for (gene, read_list) in read_dict {
        let multiple: f64 = read_list
            .iter()
            .map(|read_id| co_dict.get(&(read_id.clone(), gene)).copied().unwrap_or(0.0))
            .sum();
        *abundance.entry(gene).or_insert(0.0) += multiple;
    }

    let gene_counts: Vec<GeneCount> = abundance
        .into_iter()
        .map(|(gene_id, count)| GeneCount {
            gene_id,
            gene_length: *database.get(&gene_id).unwrap_or(&0),
            count,
        })
        .collect();
    Ok(MspCountResult {
        gene_counts,
        counted_reads,
    })
}

#[pyfunction]
fn load_fasta(path: &str) -> PyResult<HashMap<String, String>> {
    let reader = faidx::Reader::from_path(path)
        .map_err(|e| PyIOError::new_err(format!("failed to open FASTA {path}: {e}")))?;
    let names = reader
        .seq_names()
        .map_err(|e| PyIOError::new_err(format!("failed to read FASTA index for {path}: {e}")))?;

    let mut sequences = HashMap::with_capacity(names.len());
    for name in names {
        let length = reader.fetch_seq_len(&name) as usize;
        let end = length.saturating_sub(1);
        let seq = reader
            .fetch_seq_string(&name, 0, end)
            .map_err(|e| PyIOError::new_err(format!("failed to fetch sequence {name}: {e}")))?;
        sequences.insert(name, seq);
    }
    Ok(sequences)
}

#[pyfunction]
fn bed_chunks(bed_path: &str, num_chunks: usize) -> PyResult<Vec<Vec<(u32, u32, u32)>>> {
    if num_chunks == 0 {
        return Err(PyValueError::new_err("num_chunks must be > 0"));
    }

    let file = File::open(bed_path)
        .map_err(|e| PyIOError::new_err(format!("failed to open BED {bed_path}: {e}")))?;
    let reader = BufReader::new(file);

    let mut rows: Vec<(u32, u32, u32)> = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| PyIOError::new_err(format!("failed to read BED line: {e}")))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        rows.push((
            cols[0]
                .parse::<u32>()
                .map_err(|e| PyValueError::new_err(format!("invalid gene_id {}: {}", cols[0], e)))?,
            cols[1]
                .parse::<u32>()
                .map_err(|e| PyValueError::new_err(format!("invalid start {}: {}", cols[1], e)))?,
            cols[2]
                .parse::<u32>()
                .map_err(|e| PyValueError::new_err(format!("invalid end {}: {}", cols[2], e)))?,
        ));
    }

    let n = rows.len();
    let num_chunks = num_chunks.min(n.max(1));
    let base = n / num_chunks;
    let remainder = n % num_chunks;

    let mut chunks: Vec<Vec<(u32, u32, u32)>> = Vec::with_capacity(num_chunks);
    let mut start = 0;
    for i in 0..num_chunks {
        let chunk_size = base + if i < remainder { 1 } else { 0 };
        let end = start + chunk_size;
        chunks.push(rows[start..end].to_vec());
        start = end;
    }
    Ok(chunks)
}

#[pyfunction]
fn count_reads_in_gene(
    cram_path: &str,
    ref_path: &str,
    gene_name: &str,
    gene_length: u32,
    max_depth: u32,
) -> PyResult<HashMap<u32, u32>> {
    let mut reader = IndexedReader::from_path(cram_path)
        .map_err(|e| PyIOError::new_err(format!("failed to open indexed CRAM {cram_path}: {e}")))?;
    reader
        .set_reference(ref_path)
        .map_err(|e| PyIOError::new_err(format!("failed to set CRAM reference {ref_path}: {e}")))?;

    let header = reader.header().clone();
    let tid = header
        .target_names()
        .iter()
        .position(|name| bytes_to_string(name) == gene_name)
        .ok_or_else(|| PyValueError::new_err(format!("gene {gene_name} not found in CRAM header")))?;

    reader
        .fetch((tid as u32, 0, gene_length))
        .map_err(|e| PyIOError::new_err(format!("failed to fetch {gene_name}: {e}")))?;

    let mut pileup = reader.pileup();
    pileup.set_max_depth(max_depth);

    let mut counts: HashMap<u32, u32> = HashMap::new();
    for p in pileup {
        let p = p.map_err(|e| PyIOError::new_err(format!("pileup error: {e}")))?;
        let mut depth: u32 = 0;
        for alignment in p.alignments() {
            if !alignment.is_del() && !alignment.is_refskip() {
                depth += 1;
            }
        }
        counts.insert(p.pos(), depth);
    }
    Ok(counts)
}

fn iupac_code(bases: &mut Vec<char>) -> char {
    bases.sort_unstable();
    bases.dedup();
    match bases.as_slice() {
        ['A'] => 'A',
        ['C'] => 'C',
        ['G'] => 'G',
        ['T'] => 'T',
        ['A', 'G'] => 'R',
        ['A', 'C'] => 'M',
        ['A', 'T'] => 'W',
        ['C', 'G'] => 'S',
        ['C', 'T'] => 'Y',
        ['G', 'T'] => 'K',
        ['A', 'C', 'G'] => 'V',
        ['A', 'C', 'T'] => 'H',
        ['A', 'G', 'T'] => 'D',
        ['C', 'G', 'T'] => 'B',
        ['A', 'C', 'G', 'T'] => 'N',
        _ => 'N',
    }
}

#[derive(Debug)]
struct VcfVariant {
    start: u32,
    alleles: Vec<String>,
}

#[pyfunction]
fn create_consensus(
    reference_path: &str,
    vcf_path: &str,
    bed_path: &str,
    low_cov_sites: Vec<(u32, u32, u32)>,
    gene_ignore: Vec<(u32, u32)>,
    min_frequency: f64,
    gap_char: char,
) -> PyResult<Vec<(u32, String)>> {
    let ref_reader = faidx::Reader::from_path(reference_path)
        .map_err(|e| PyIOError::new_err(format!("failed to open FASTA {reference_path}: {e}")))?;
    let ref_names = ref_reader
        .seq_names()
        .map_err(|e| PyIOError::new_err(format!("failed to read FASTA index for {reference_path}: {e}")))?;
    let mut references: HashMap<u32, String> = HashMap::with_capacity(ref_names.len());
    for name in ref_names {
        let gene_id = name.parse::<u32>().map_err(|_| {
            PyValueError::new_err(format!("reference name {name} is not a numeric gene id"))
        })?;
        let length = ref_reader.fetch_seq_len(&name) as usize;
        let end = length.saturating_sub(1);
        let seq = ref_reader
            .fetch_seq_string(&name, 0, end)
            .map_err(|e| PyIOError::new_err(format!("failed to fetch reference {name}: {e}")))?;
        references.insert(gene_id, seq);
    }

    let bed_file = File::open(bed_path)
        .map_err(|e| PyIOError::new_err(format!("failed to open BED {bed_path}: {e}")))?;
    let bed_reader = BufReader::new(bed_file);
    let mut gene_ids: Vec<u32> = Vec::new();
    for line in bed_reader.lines() {
        let line = line.map_err(|e| PyIOError::new_err(format!("failed to read BED line: {e}")))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(gene_str) = line.split('\t').next() {
            gene_ids.push(
                gene_str
                    .parse::<u32>()
                    .map_err(|e| PyValueError::new_err(format!("invalid gene_id {gene_str}: {e}")))?,
            );
        }
    }
    gene_ids.sort_unstable();
    gene_ids.dedup();

    let mut low_cov_map: HashMap<u32, Vec<(u32, u32)>> = HashMap::new();
    for (gene_id, start, end) in low_cov_sites {
        low_cov_map.entry(gene_id).or_default().push((start, end));
    }
    let mut ignore_map: HashMap<u32, u32> = HashMap::new();
    for (gene_id, length) in gene_ignore {
        ignore_map.insert(gene_id, length);
    }

    let mut variants_by_gene: HashMap<u32, Vec<VcfVariant>> = HashMap::new();
    {
        use rust_htslib::bcf::Read;
        let mut vcf_reader = rust_htslib::bcf::Reader::from_path(vcf_path)
            .map_err(|e| PyIOError::new_err(format!("failed to open VCF {vcf_path}: {e}")))?;
        let header = vcf_reader.header().clone();
        for record in vcf_reader.records() {
            let record = record.map_err(|e| PyIOError::new_err(format!("VCF read error: {e}")))?;
            let rid = record.rid().ok_or_else(|| {
                PyValueError::new_err("VCF record missing chromosome id")
            })?;
            let chrom = String::from_utf8_lossy(header.rid2name(rid).map_err(|e| {
                PyValueError::new_err(format!("VCF header rid lookup failed: {e}"))
            })?)
            .to_string();
            let gene_id = chrom.parse::<u32>().map_err(|_| {
                PyValueError::new_err(format!("VCF chromosome {chrom} is not a numeric gene id"))
            })?;

            let start = record.pos() as u32;

            let alleles: Vec<String> = record
                .alleles()
                .into_iter()
                .map(|a| String::from_utf8_lossy(a).to_string().to_ascii_uppercase())
                .collect();
            if alleles.is_empty() {
                continue;
            }

            let ro = record
                .info(b"RO")
                .integer()
                .map_err(|e| PyValueError::new_err(format!("VCF INFO/RO read error: {e}")))?
                .and_then(|v| v.first().copied())
                .unwrap_or(0) as f64;
            let ao_sum = record
                .info(b"AO")
                .integer()
                .map_err(|e| PyValueError::new_err(format!("VCF INFO/AO read error: {e}")))?
                .map(|v| v.iter().map(|&x| x as f64).sum::<f64>())
                .unwrap_or(0.0);
            let denominator = ro + ao_sum;
            let ref_frequency = if denominator > 0.0 { ro / denominator } else { 0.0 };

            let keep_alleles: Vec<String> = if ref_frequency >= min_frequency {
                alleles
            } else {
                alleles.into_iter().skip(1).collect()
            };
            if keep_alleles.is_empty() {
                continue;
            }

            variants_by_gene.entry(gene_id).or_default().push(VcfVariant {
                start,
                alleles: keep_alleles,
            });
        }
    }

    let mut results: Vec<(u32, String)> = Vec::with_capacity(gene_ids.len());
    for gene_id in gene_ids {
        let consensus = if let Some(&length) = ignore_map.get(&gene_id) {
            gap_char.to_string().repeat(length as usize)
        } else {
            let ref_seq = references.get(&gene_id).ok_or_else(|| {
                PyValueError::new_err(format!("reference sequence for gene {gene_id} not found"))
            })?;
            let mut chars: Vec<char> = ref_seq.chars().collect();
            if let Some(variants) = variants_by_gene.get(&gene_id) {
                for variant in variants {
                    let max_len = variant
                        .alleles
                        .iter()
                        .map(|a| a.len())
                        .max()
                        .unwrap_or(1);
                    for i in 0..max_len {
                        let mut bases: Vec<char> = variant
                            .alleles
                            .iter()
                            .filter_map(|a| a.chars().nth(i))
                            .collect();
                        let pos = variant.start + i as u32;
                        if pos < chars.len() as u32 && !bases.is_empty() {
                            chars[pos as usize] = iupac_code(&mut bases);
                        }
                    }
                }
            }
            if let Some(ranges) = low_cov_map.get(&gene_id) {
                for &(start, end) in ranges {
                    let start = start as usize;
                    let end = (end as usize).min(chars.len());
                    for c in chars.iter_mut().take(end).skip(start) {
                        *c = gap_char;
                    }
                }
            }
            chars.into_iter().collect()
        };
        results.push((gene_id, consensus));
    }
    Ok(results)
}

#[pymodule]
fn meteor_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(count_cram_records, m)?)?;
    m.add_function(wrap_pyfunction!(sum_cram_cigar_lengths, m)?)?;
    m.add_function(wrap_pyfunction!(cram_records, m)?)?;
    m.add_function(wrap_pyfunction!(count_msp, m)?)?;
    m.add_function(wrap_pyfunction!(load_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(bed_chunks, m)?)?;
    m.add_function(wrap_pyfunction!(count_reads_in_gene, m)?)?;
    m.add_function(wrap_pyfunction!(create_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(freebayes::call_variants_parallel, m)?)?;
    m.add_function(wrap_pyfunction!(vcf::write_vcf_text, m)?)?;
    m.add_function(wrap_pyfunction!(vcf::bgzip_file, m)?)?;
    m.add_function(wrap_pyfunction!(vcf::write_bcf, m)?)?;
    m.add_class::<CramRecord>()?;
    m.add_class::<GeneCount>()?;
    m.add_class::<MspCountResult>()?;
    m.add_class::<freebayes::FreebayesOptions>()?;
    m.add("FreebayesError", m.py().get_type::<freebayes::FreebayesError>())?;
    m.add_class::<vcf::VcfRecord>()?;
    m.add("VcfError", m.py().get_type::<vcf::VcfError>())?;
    Ok(())
}
