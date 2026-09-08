use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;

pyo3::create_exception!(meteor_core, VcfError, pyo3::exceptions::PyException);

#[pyclass]
#[derive(Clone, Debug)]
pub struct VcfRecord {
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub pos: u32,
    #[pyo3(get)]
    pub id: String,
    #[pyo3(get)]
    pub ref_allele: String,
    #[pyo3(get)]
    pub alt_alleles: Vec<String>,
    #[pyo3(get)]
    pub qual: f32,
    #[pyo3(get)]
    pub filter: String,
    #[pyo3(get)]
    pub info: HashMap<String, String>,
}

#[pymethods]
impl VcfRecord {
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        chrom: String,
        pos: u32,
        id: String,
        ref_allele: String,
        alt_alleles: Vec<String>,
        qual: f32,
        filter: String,
        info: HashMap<String, String>,
    ) -> PyResult<Self> {
        if chrom.is_empty() {
            return Err(VcfError::new_err("chrom must not be empty"));
        }
        if pos == 0 {
            return Err(VcfError::new_err(
                "VCF positions are 1-based and must be > 0",
            ));
        }
        if ref_allele.is_empty() {
            return Err(VcfError::new_err("ref_allele must not be empty"));
        }
        if !ref_allele
            .chars()
            .all(|c| matches!(c, 'A' | 'C' | 'G' | 'T' | 'N' | 'a' | 'c' | 'g' | 't' | 'n'))
        {
            return Err(VcfError::new_err(format!(
                "ref_allele contains invalid bases: {ref_allele}"
            )));
        }
        for alt in &alt_alleles {
            if alt != "."
                && !alt.chars().all(|c| {
                    matches!(
                        c,
                        'A' | 'C' | 'G' | 'T' | 'N' | 'a' | 'c' | 'g' | 't' | 'n' | '*'
                    )
                })
            {
                return Err(VcfError::new_err(format!(
                    "alt_allele contains invalid bases: {alt}"
                )));
            }
        }
        Ok(Self {
            chrom,
            pos,
            id,
            ref_allele,
            alt_alleles,
            qual,
            filter,
            info,
        })
    }
}

fn contig_lengths(records: &[VcfRecord]) -> HashMap<String, u32> {
    let mut lengths: HashMap<String, u32> = HashMap::new();
    for rec in records {
        let max_allele_len = rec
            .alt_alleles
            .iter()
            .chain(std::iter::once(&rec.ref_allele))
            .map(String::len)
            .max()
            .unwrap_or(1) as u32;
        let end = rec.pos + max_allele_len - 1;
        let entry = lengths.entry(rec.chrom.clone()).or_insert(0);
        if end > *entry {
            *entry = end;
        }
    }
    lengths
}

fn format_info(info: &HashMap<String, String>) -> String {
    if info.is_empty() {
        return ".".to_string();
    }
    let mut pairs: Vec<String> = Vec::with_capacity(info.len());
    for (k, v) in info {
        if v.is_empty() {
            pairs.push(k.clone());
        } else {
            pairs.push(format!("{k}={v}"));
        }
    }
    pairs.sort_unstable();
    pairs.join(";")
}

#[pyfunction]
pub fn write_vcf_text(
    records: Vec<VcfRecord>,
    output_path: &str,
    sample_name: &str,
) -> PyResult<()> {
    let mut file = File::create(output_path)
        .map_err(|e| PyIOError::new_err(format!("failed to create VCF {output_path}: {e}")))?;

    writeln!(file, "##fileformat=VCFv4.2")
        .map_err(|e| PyIOError::new_err(format!("failed to write VCF header: {e}")))?;
    writeln!(file, "##source=meteor-core")
        .map_err(|e| PyIOError::new_err(format!("failed to write VCF header: {e}")))?;

    let mut info_keys: HashSet<String> = HashSet::new();
    for rec in &records {
        for key in rec.info.keys() {
            info_keys.insert(key.clone());
        }
    }
    for key in info_keys {
        writeln!(
            file,
            "##INFO=<ID={key},Number=.,Type=String,Description=\"{key}\">"
        )
        .map_err(|e| PyIOError::new_err(format!("failed to write VCF header: {e}")))?;
    }
    writeln!(
        file,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )
    .map_err(|e| PyIOError::new_err(format!("failed to write VCF header: {e}")))?;

    for (chrom, length) in contig_lengths(&records) {
        writeln!(file, "##contig=<ID={chrom},length={length}>")
            .map_err(|e| PyIOError::new_err(format!("failed to write VCF header: {e}")))?;
    }

    writeln!(
        file,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"
    )
    .map_err(|e| PyIOError::new_err(format!("failed to write VCF header: {e}")))?;

    for rec in records {
        let alt = if rec.alt_alleles.is_empty() {
            ".".to_string()
        } else {
            rec.alt_alleles.join(",")
        };
        let qual = if rec.qual < 0.0 {
            ".".to_string()
        } else {
            rec.qual.to_string()
        };
        let filter = if rec.filter.is_empty() {
            ".".to_string()
        } else {
            rec.filter
        };
        let info = format_info(&rec.info);
        writeln!(
            file,
            "{chrom}\t{pos}\t{id}\t{ref_allele}\t{alt}\t{qual}\t{filter}\t{info}\tGT\t./.",
            chrom = rec.chrom,
            pos = rec.pos,
            id = if rec.id.is_empty() {
                ".".to_string()
            } else {
                rec.id
            },
            ref_allele = rec.ref_allele,
        )
        .map_err(|e| PyIOError::new_err(format!("failed to write VCF record: {e}")))?;
    }

    Ok(())
}

#[pyfunction]
pub fn bgzip_file(input_path: &str, output_path: &str) -> PyResult<()> {
    let reader = File::open(input_path)
        .map_err(|e| PyIOError::new_err(format!("failed to open input {input_path}: {e}")))?;
    let reader = BufReader::new(reader);
    let mut writer = rust_htslib::bgzf::Writer::from_path(output_path).map_err(|e| {
        PyIOError::new_err(format!("failed to create bgzip output {output_path}: {e}"))
    })?;

    for line in reader.lines() {
        let line =
            line.map_err(|e| PyIOError::new_err(format!("failed to read input line: {e}")))?;
        writer
            .write_all(line.as_bytes())
            .map_err(|e| PyIOError::new_err(format!("failed to write bgzip data: {e}")))?;
        writer
            .write_all(b"\n")
            .map_err(|e| PyIOError::new_err(format!("failed to write bgzip data: {e}")))?;
    }
    writer
        .flush()
        .map_err(|e| PyIOError::new_err(format!("failed to flush bgzip writer: {e}")))?;
    Ok(())
}

#[pyfunction]
pub fn write_bcf(records: Vec<VcfRecord>, output_path: &str, sample_name: &str) -> PyResult<()> {
    let temp_dir = tempfile::Builder::new()
        .prefix("meteor_bcf_")
        .tempdir()
        .map_err(|e| PyIOError::new_err(format!("failed to create temp dir: {e}")))?;
    let temp_vcf = temp_dir.path().join("temp.vcf");
    write_vcf_text(
        records,
        temp_vcf.to_str().expect("temp path is valid UTF-8"),
        sample_name,
    )?;

    use rust_htslib::bcf::Read;
    let mut reader = rust_htslib::bcf::Reader::from_path(&temp_vcf).map_err(|e| {
        PyIOError::new_err(format!("failed to open temp VCF for BCF conversion: {e}"))
    })?;
    let header = rust_htslib::bcf::Header::from_template(reader.header());
    let mut writer = rust_htslib::bcf::Writer::from_path(
        output_path,
        &header,
        false,
        rust_htslib::bcf::Format::Bcf,
    )
    .map_err(|e| PyIOError::new_err(format!("failed to create BCF writer {output_path}: {e}")))?;

    for record in reader.records() {
        let record = record
            .map_err(|e| PyIOError::new_err(format!("failed to read temp VCF record: {e}")))?;
        writer
            .write(&record)
            .map_err(|e| PyIOError::new_err(format!("failed to write BCF record: {e}")))?;
    }
    Ok(())
}
