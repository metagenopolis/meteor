use std::collections::HashMap;

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{Read, Reader, Record};

use crate::threads;

#[pyclass]
#[derive(Clone)]
pub struct GeneDepth {
    #[pyo3(get)]
    pub gene: String,
    #[pyo3(get)]
    pub depths: Vec<u32>,
}

struct GeneInterval {
    start: u64,
    end: u64,
    index: usize,
}

pub(crate) fn open_cram(cram_path: &str, ref_path: Option<&str>) -> PyResult<Reader> {
    let mut reader = Reader::from_path(cram_path)
        .map_err(|e| PyIOError::new_err(format!("failed to open CRAM {cram_path}: {e}")))?;
    if let Some(path) = ref_path {
        reader
            .set_reference(path)
            .map_err(|e| PyIOError::new_err(format!("failed to set CRAM reference {path}: {e}")))?;
    }
    reader
        .set_threads(threads::resolve_thread_count())
        .map_err(|e| PyIOError::new_err(format!("failed to set CRAM thread pool: {e}")))?;
    Ok(reader)
}

fn gene_tid(header: &rust_htslib::bam::HeaderView, gene_name: &str) -> Option<usize> {
    header
        .target_names()
        .iter()
        .position(|name| crate::bytes_to_string(name) == gene_name)
}

fn skip_read(record: &Record) -> bool {
    record.is_unmapped()
        || record.is_secondary()
        || record.is_duplicate()
        || record.is_quality_check_failed()
}

fn add_coverage(
    intervals: &[GeneInterval],
    block_start: u64,
    block_end: u64,
    depths: &mut [Vec<u32>],
    max_depth: u32,
    active_start: &mut usize,
) {
    if block_start >= block_end {
        return;
    }
    while *active_start < intervals.len() && intervals[*active_start].end <= block_start {
        *active_start += 1;
    }
    let end_idx = intervals.partition_point(|iv| iv.start < block_end);
    for iv in &intervals[*active_start..end_idx] {
        if iv.end <= block_start {
            continue;
        }
        let overlap_start = block_start.max(iv.start);
        let overlap_end = block_end.min(iv.end);
        if overlap_start >= overlap_end {
            continue;
        }
        let start_idx = (overlap_start - iv.start) as usize;
        let end_idx_rel = (overlap_end - iv.start) as usize;
        let iv_depths = &mut depths[iv.index];
        for depth in &mut iv_depths[start_idx..end_idx_rel] {
            if *depth < max_depth {
                *depth += 1;
            }
        }
    }
}

fn walk_cigar(
    record: &Record,
    intervals: &[GeneInterval],
    depths: &mut [Vec<u32>],
    max_depth: u32,
) {
    let mut ref_pos = record.pos();
    let mut active_start = 0usize;
    for cigar in record.cigar().iter() {
        match cigar {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = i64::from(*len);
                add_coverage(
                    intervals,
                    ref_pos as u64,
                    (ref_pos + len_i64) as u64,
                    depths,
                    max_depth,
                    &mut active_start,
                );
                ref_pos += len_i64;
            }
            Cigar::Ins(_) => {}
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += i64::from(*len);
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
}

#[pyfunction]
pub fn depth_per_gene(
    cram_path: &str,
    ref_path: &str,
    genes: Vec<(String, u64, u64)>,
    max_depth: u32,
) -> PyResult<Vec<GeneDepth>> {
    if genes.is_empty() {
        return Ok(Vec::new());
    }

    let mut reader = open_cram(cram_path, Some(ref_path))?;
    let header = reader.header().clone();

    let mut depths: Vec<Vec<u32>> = Vec::with_capacity(genes.len());
    let mut gene_names: Vec<String> = Vec::with_capacity(genes.len());
    let mut by_tid: HashMap<i32, Vec<GeneInterval>> = HashMap::new();

    for (index, (name, start, end)) in genes.into_iter().enumerate() {
        if start >= end {
            return Err(PyValueError::new_err(format!(
                "invalid interval for {name}: start {start} >= end {end}"
            )));
        }
        let tid = gene_tid(&header, &name)
            .ok_or_else(|| PyValueError::new_err(format!("gene {name} not found in CRAM header")))?
            as i32;
        let contig_len = header.target_len(tid as u32).unwrap_or(0);
        let end = end.min(contig_len);
        let len = (end - start) as usize;
        depths.push(vec![0; len]);
        gene_names.push(name);
        by_tid
            .entry(tid)
            .or_default()
            .push(GeneInterval { start, end, index });
    }

    for intervals in by_tid.values_mut() {
        intervals.sort_by_key(|iv| iv.start);
    }

    for record in reader.records() {
        let record = record.map_err(|e| PyIOError::new_err(format!("CRAM read error: {e}")))?;
        if skip_read(&record) {
            continue;
        }
        let tid = record.tid();
        if tid < 0 {
            continue;
        }
        let intervals = match by_tid.get(&tid) {
            Some(ivs) => ivs,
            None => continue,
        };
        walk_cigar(&record, intervals, &mut depths, max_depth);
    }

    let mut results = Vec::with_capacity(gene_names.len());
    for (gene, gene_depths) in gene_names.into_iter().zip(depths) {
        results.push(GeneDepth {
            gene,
            depths: gene_depths,
        });
    }
    Ok(results)
}
