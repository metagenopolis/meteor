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

#[pymodule]
fn meteor_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(count_cram_records, m)?)?;
    m.add_function(wrap_pyfunction!(sum_cram_cigar_lengths, m)?)?;
    m.add_function(wrap_pyfunction!(cram_records, m)?)?;
    m.add_class::<CramRecord>()?;
    Ok(())
}
