use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use rust_htslib::bam::{Read, Reader, Record};

use crate::{aligned_nucleotides, bytes_to_string, extract_nm};

/// One per-gene aggregate returned by ``count_msp_aggregates``.
///
/// ``reads`` is a pre-joined newline-delimited string of the read ids that
/// contributed to this gene. Keeping it as one string collapses the Rust/Python
/// boundary to roughly one crossing per gene instead of one per read-gene pair.
#[pyclass]
#[derive(Clone)]
pub struct AggregateRow {
    #[pyo3(get)]
    pub msp: String,
    #[pyo3(get)]
    pub gene: String,
    #[pyo3(get)]
    pub count: f64,
    #[pyo3(get)]
    pub reads: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountingType {
    SmartShared,
    Unique,
    Total,
}

impl CountingType {
    pub fn from_str(counting_type: &str) -> PyResult<Self> {
        match counting_type {
            "smart_shared" => Ok(Self::SmartShared),
            "unique" => Ok(Self::Unique),
            "total" => Ok(Self::Total),
            _ => Err(PyErr::new::<PyValueError, _>(format!(
                "{counting_type} is not a valid counting type"
            ))),
        }
    }
}

/// Shared result type produced by the counting core.
pub struct CountCoreResult {
    pub database: BTreeMap<i32, i32>,
    pub abundance: BTreeMap<i32, f64>,
    pub gene_reads: BTreeMap<i32, Vec<String>>,
    pub counted_reads: usize,
}

/// Open a CRAM file, inferring the reference FASTA from the MSP map directory.
fn open_cram_with_msp_map(cram_path: &str, msp_map_path: &str) -> PyResult<Reader> {
    let mut reader = Reader::from_path(cram_path)
        .map_err(|e| PyIOError::new_err(format!("failed to open CRAM {cram_path}: {e}")))?;
    if let Some(ref_path) = find_reference_path(msp_map_path) {
        reader
            .set_reference(
                ref_path
                    .to_str()
                    .ok_or_else(|| PyValueError::new_err("reference path is not valid UTF-8"))?,
            )
            .map_err(|e| {
                PyIOError::new_err(format!(
                    "failed to set CRAM reference {}: {e}",
                    ref_path.display()
                ))
            })?;
    }
    Ok(reader)
}

/// Look for a reference FASTA next to the MSP map file.
fn find_reference_path(msp_map_path: &str) -> Option<PathBuf> {
    let base = Path::new(msp_map_path).parent()?;
    for name in [
        "reference.fa",
        "reference.fa.gz",
        "reference.fasta",
        "reference.fasta.gz",
    ] {
        let candidate = base.join(name);
        if candidate.is_file() {
            return Some(candidate);
        }
    }
    None
}

/// Parse the MSP map TSV into a gene_id -> msp_name lookup.
fn load_msp_map(path: &str) -> PyResult<HashMap<i32, String>> {
    let file = File::open(path)
        .map_err(|e| PyIOError::new_err(format!("failed to open MSP map {path}: {e}")))?;
    let reader = BufReader::new(file);
    let mut mapping = HashMap::new();

    for (idx, line) in reader.lines().enumerate() {
        let line = line
            .map_err(|e| PyIOError::new_err(format!("failed to read MSP map line {idx}: {e}")))?;
        if idx == 0 || line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            continue;
        }
        let gene_id = cols[1]
            .parse::<i32>()
            .map_err(|e| PyValueError::new_err(format!("invalid gene_id {}: {}", cols[1], e)))?;
        mapping.insert(gene_id, cols[0].to_string());
    }
    Ok(mapping)
}

/// Shared CRAM iteration and counting core used by ``count_msp`` and
/// ``count_msp_aggregates``.
///
/// The algorithm is intentionally identical to ``meteor.counter.Counter`` so that
/// the Rust and Python paths produce byte-identical TSV files.
pub fn count_msp_core(
    reader: &mut Reader,
    identity_threshold: f64,
    counting_type: CountingType,
) -> PyResult<CountCoreResult> {
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
    let mut gene_reads: HashMap<i32, Vec<String>> = HashMap::new();

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
                reads.insert(query_name.clone(), (score, vec![gene_id]));
            }
        }
        gene_reads.entry(gene_id).or_default().push(query_name);
    }

    let counted_reads = reads.len();

    let mut gene_reads: BTreeMap<i32, Vec<String>> = gene_reads
        .into_iter()
        .map(|(gene, mut names)| {
            names.sort_unstable();
            names.dedup();
            (gene, names)
        })
        .collect();

    if counting_type == CountingType::Total {
        let mut abundance: BTreeMap<i32, f64> = database.keys().map(|&g| (g, 0.0)).collect();
        for genes in reads.values() {
            for gene in &genes.1 {
                *abundance.entry(*gene).or_insert(0.0) += 1.0;
            }
        }
        for &gene in database.keys() {
            gene_reads.entry(gene).or_default();
        }
        return Ok(CountCoreResult {
            database,
            abundance,
            gene_reads,
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

    if counting_type == CountingType::Unique {
        for &gene in database.keys() {
            gene_reads.entry(gene).or_default();
        }
        return Ok(CountCoreResult {
            database,
            abundance: unique_on_gene,
            gene_reads,
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

        let duplicated_genes: HashSet<i32> = genes
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
            .map(|read_id| {
                co_dict
                    .get(&(read_id.clone(), gene))
                    .copied()
                    .unwrap_or(0.0)
            })
            .sum();
        *abundance.entry(gene).or_insert(0.0) += multiple;
    }

    for &gene in database.keys() {
        gene_reads.entry(gene).or_default();
    }

    Ok(CountCoreResult {
        database,
        abundance,
        gene_reads,
        counted_reads,
    })
}

/// Count reads per MSP/gene and return aggregate rows.
///
/// This is the aggregates-only replacement for ``count_msp``: instead of
/// pushing per-read strings across the PyO3 boundary, each gene crosses once
/// with its pre-joined read list.
#[pyfunction]
pub fn count_msp_aggregates(
    cram_path: &str,
    msp_map_path: &str,
    identity_threshold: f64,
    counting_type: &str,
) -> PyResult<Vec<AggregateRow>> {
    let counting_type = CountingType::from_str(counting_type)?;
    let mut reader = open_cram_with_msp_map(cram_path, msp_map_path)?;
    let core = count_msp_core(&mut reader, identity_threshold, counting_type)?;
    let msp_map = load_msp_map(msp_map_path)?;

    let mut rows: Vec<AggregateRow> = Vec::with_capacity(core.database.len());
    for (gene_id, _) in core.database {
        let count = core.abundance.get(&gene_id).copied().unwrap_or(0.0);
        let reads = core
            .gene_reads
            .get(&gene_id)
            .map(|names| names.join("\n"))
            .unwrap_or_default();
        let msp = msp_map.get(&gene_id).cloned().unwrap_or_default();
        rows.push(AggregateRow {
            msp,
            gene: gene_id.to_string(),
            count,
            reads,
        });
    }
    Ok(rows)
}
