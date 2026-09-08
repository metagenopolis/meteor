use pyo3::prelude::*;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{Read, Reader, Record};

use crate::{bytes_to_string, extract_nm, open_cram};

#[pyclass]
#[derive(Clone)]
pub struct CramRecord {
    #[pyo3(get)]
    pub query_name: String,
    #[pyo3(get)]
    pub reference_name: String,
    #[pyo3(get)]
    pub cigar_tuples: Vec<(u32, u32)>,
    #[pyo3(get)]
    pub nm: Option<u32>,
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

fn build_cram_record(record: &Record, target_names: &[String]) -> CramRecord {
    let query_name = bytes_to_string(record.qname());
    let reference_name = if record.tid() >= 0 && (record.tid() as usize) < target_names.len() {
        target_names[record.tid() as usize].clone()
    } else {
        "*".to_string()
    };
    CramRecord {
        query_name,
        reference_name,
        cigar_tuples: cigar_tuples(record),
        nm: extract_nm(record),
    }
}

#[pyfunction]
pub fn count_cram_records(cram_path: &str, ref_path: &str) -> PyResult<usize> {
    let mut reader = open_cram(cram_path, ref_path)?;
    let mut record = Record::new();
    let mut count: usize = 0;
    while let Some(result) = reader.read(&mut record) {
        result
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("CRAM read error: {e}")))?;
        count += 1;
    }
    Ok(count)
}

#[pyfunction]
pub fn sum_cram_cigar_lengths(cram_path: &str, ref_path: &str) -> PyResult<(u64, u64)> {
    let mut reader = open_cram(cram_path, ref_path)?;
    let mut record = Record::new();
    let mut total_len: u64 = 0;
    let mut total_ops: u64 = 0;
    while let Some(result) = reader.read(&mut record) {
        result
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("CRAM read error: {e}")))?;
        for cigar in record.cigar().iter() {
            total_len += u64::from(cigar.len());
            total_ops += 1;
        }
    }
    Ok((total_len, total_ops))
}

#[pyfunction]
pub fn cram_records(cram_path: &str, ref_path: &str) -> PyResult<Vec<CramRecord>> {
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
        result
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("CRAM read error: {e}")))?;
        records.push(build_cram_record(&record, &target_names));
    }
    Ok(records)
}

#[pyclass(unsendable)]
pub struct CramRecordStream {
    reader: Option<Reader>,
    target_names: Vec<String>,
    buffered: Vec<CramRecord>,
    buffered_index: usize,
}

#[pymethods]
impl CramRecordStream {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<CramRecord>> {
        if slf.buffered_index < slf.buffered.len() {
            let rec = slf.buffered[slf.buffered_index].clone();
            slf.buffered_index += 1;
            return Ok(Some(rec));
        }
        let reader = slf
            .reader
            .as_mut()
            .ok_or_else(|| pyo3::exceptions::PyIOError::new_err("CRAM reader is not available"))?;
        let mut record = Record::new();
        match reader.read(&mut record) {
            Some(result) => {
                result.map_err(|e| {
                    pyo3::exceptions::PyIOError::new_err(format!("CRAM read error: {e}"))
                })?;
                Ok(Some(build_cram_record(&record, &slf.target_names)))
            }
            None => Ok(None),
        }
    }
}

#[pyfunction]
pub fn stream_cram_records(
    py: Python<'_>,
    cram_path: &str,
    ref_path: &str,
) -> PyResult<CramRecordStream> {
    match open_cram(cram_path, ref_path) {
        Ok(reader) => {
            let header = reader.header().clone();
            let target_names: Vec<String> = header
                .target_names()
                .iter()
                .map(|name| bytes_to_string(name))
                .collect();
            Ok(CramRecordStream {
                reader: Some(reader),
                target_names,
                buffered: Vec::new(),
                buffered_index: 0,
            })
        }
        Err(e) => {
            let warning = py
                .import("logging")
                .and_then(|m| m.getattr("warning"))
                .map_err(|_| {
                    pyo3::exceptions::PyRuntimeError::new_err("failed to import logging")
                })?;
            let _ = warning.call1((format!(
                "Streaming CRAM read failed ({e}); falling back to buffered mode for {cram_path}"
            ),));
            let buffered = cram_records(cram_path, ref_path)?;
            Ok(CramRecordStream {
                reader: None,
                target_names: Vec::new(),
                buffered,
                buffered_index: 0,
            })
        }
    }
}
