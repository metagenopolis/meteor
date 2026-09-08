use std::collections::{BTreeMap, HashMap};

use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{Read, Reader, Record};

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
        Aux::U32(v) => Some(v as u32),
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
        let name = bytes_to_string(&header.target_names()[tid as usize]);
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
            bytes_to_string(&header.target_names()[record.tid() as usize])
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

#[pymodule]
fn meteor_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(count_cram_records, m)?)?;
    m.add_function(wrap_pyfunction!(sum_cram_cigar_lengths, m)?)?;
    m.add_function(wrap_pyfunction!(cram_records, m)?)?;
    m.add_function(wrap_pyfunction!(count_msp, m)?)?;
    m.add_class::<CramRecord>()?;
    m.add_class::<GeneCount>()?;
    m.add_class::<MspCountResult>()?;
    Ok(())
}
