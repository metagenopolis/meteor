use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::process::{Command, Stdio};
use std::sync::mpsc::{channel, Sender};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use wait_timeout::ChildExt;

pyo3::create_exception!(meteor_core, FreebayesError, pyo3::exceptions::PyException);

#[pyclass]
#[derive(Clone)]
pub struct FreebayesOptions {
    #[pyo3(get)]
    pub min_snp_depth: u32,
    #[pyo3(get)]
    pub min_frequency: f64,
    #[pyo3(get)]
    pub ploidy: u32,
}

#[pymethods]
impl FreebayesOptions {
    #[new]
    fn new(min_snp_depth: u32, min_frequency: f64, ploidy: u32) -> Self {
        Self {
            min_snp_depth,
            min_frequency,
            ploidy,
        }
    }
}

struct ChunkTask {
    index: usize,
    bed_path: std::path::PathBuf,
    stdout_path: std::path::PathBuf,
    stderr_path: std::path::PathBuf,
}

enum ChunkOutcome {
    Success { index: usize, stdout_path: std::path::PathBuf },
    Failure {
        index: usize,
        exit_code: Option<i32>,
        stderr: String,
    },
}

fn run_freebayes_chunk(
    freebayes_path: &str,
    cram_path: &str,
    fasta_path: &str,
    options: &FreebayesOptions,
    task: &ChunkTask,
    timeout: Duration,
) -> ChunkOutcome {
    let stdout_file = match File::create(&task.stdout_path) {
        Ok(f) => f,
        Err(e) => {
            return ChunkOutcome::Failure {
                index: task.index,
                exit_code: None,
                stderr: format!("failed to create stdout file: {e}"),
            };
        }
    };
    let stderr_file = match File::create(&task.stderr_path) {
        Ok(f) => f,
        Err(e) => {
            return ChunkOutcome::Failure {
                index: task.index,
                exit_code: None,
                stderr: format!("failed to create stderr file: {e}"),
            };
        }
    };

    let mut child = match Command::new(freebayes_path)
        .arg("-i")
        .arg("-u")
        .arg("--pooled-continuous")
        .arg("--haplotype-length")
        .arg("0")
        .arg("--min-alternate-count")
        .arg("1")
        .arg("--min-coverage")
        .arg(options.min_snp_depth.to_string())
        .arg("--min-alternate-fraction")
        .arg(options.min_frequency.to_string())
        .arg("--min-mapping-quality")
        .arg("0")
        .arg("--use-duplicate-reads")
        .arg("-t")
        .arg(&task.bed_path)
        .arg("-p")
        .arg(options.ploidy.to_string())
        .arg("-f")
        .arg(fasta_path)
        .arg("-b")
        .arg(cram_path)
        .stdout(Stdio::from(stdout_file))
        .stderr(Stdio::from(stderr_file))
        .spawn()
    {
        Ok(c) => c,
        Err(e) => {
            return ChunkOutcome::Failure {
                index: task.index,
                exit_code: None,
                stderr: format!("failed to spawn freebayes: {e}"),
            };
        }
    };

    let status = match child.wait_timeout(timeout) {
        Ok(Some(s)) => s,
        Ok(None) => {
            let _ = child.kill();
            let _ = child.wait();
            return ChunkOutcome::Failure {
                index: task.index,
                exit_code: None,
                stderr: format!("freebayes timed out after {timeout:?}"),
            };
        }
        Err(e) => {
            return ChunkOutcome::Failure {
                index: task.index,
                exit_code: None,
                stderr: format!("error waiting for freebayes: {e}"),
            };
        }
    };

    let mut stderr = String::new();
    if let Ok(mut f) = File::open(&task.stderr_path) {
        let _ = f.read_to_string(&mut stderr);
    }

    if status.success() {
        ChunkOutcome::Success {
            index: task.index,
            stdout_path: task.stdout_path.clone(),
        }
    } else {
        ChunkOutcome::Failure {
            index: task.index,
            exit_code: status.code(),
            stderr,
        }
    }
}

fn read_bed_lines(bed_path: &str) -> PyResult<Vec<String>> {
    let file = File::open(bed_path)
        .map_err(|e| PyIOError::new_err(format!("failed to open BED {bed_path}: {e}")))?;
    let reader = BufReader::new(file);
    let mut lines = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| PyIOError::new_err(format!("failed to read BED line: {e}")))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        lines.push(trimmed.to_string());
    }
    Ok(lines)
}

fn merge_vcf_chunks(stdout_paths: &[std::path::PathBuf]) -> PyResult<String> {
    let mut merged = String::new();
    let mut header_emitted = false;
    for path in stdout_paths {
        let file = File::open(path)
            .map_err(|e| PyIOError::new_err(format!("failed to open chunk VCF {path:?}: {e}")))?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line =
                line.map_err(|e| PyIOError::new_err(format!("failed to read chunk VCF line: {e}")))?;
            if line.starts_with('#') {
                if !header_emitted {
                    merged.push_str(&line);
                    merged.push('\n');
                }
            } else {
                header_emitted = true;
                merged.push_str(&line);
                merged.push('\n');
            }
        }
        header_emitted = true;
    }
    Ok(merged)
}

#[pyfunction]
pub fn call_variants_parallel(
    cram_path: &str,
    fasta_path: &str,
    bed_path: &str,
    freebayes_path: &str,
    options: &FreebayesOptions,
    n_threads: usize,
    output_path: Option<&str>,
) -> PyResult<String> {
    let bed_lines = read_bed_lines(bed_path)?;
    if bed_lines.is_empty() {
        return Err(FreebayesError::new_err("BED file contains no regions"));
    }

    let n_threads = n_threads.max(1);
    let num_chunks = n_threads.min(bed_lines.len());
    let base = bed_lines.len() / num_chunks;
    let remainder = bed_lines.len() % num_chunks;

    let temp_dir = tempfile::Builder::new()
        .prefix("meteor_freebayes_")
        .tempdir()
        .map_err(|e| PyIOError::new_err(format!("failed to create temp dir: {e}")))?;

    let mut chunks: Vec<ChunkTask> = Vec::with_capacity(num_chunks);
    let mut start = 0;
    for i in 0..num_chunks {
        let chunk_size = base + if i < remainder { 1 } else { 0 };
        let end = start + chunk_size;
        let bed_chunk_path = temp_dir.path().join(format!("chunk_{i}.bed"));
        let stdout_path = temp_dir.path().join(format!("chunk_{i}.vcf"));
        let stderr_path = temp_dir.path().join(format!("chunk_{i}.err"));

        let mut bed_file = File::create(&bed_chunk_path)
            .map_err(|e| PyIOError::new_err(format!("failed to create chunk BED: {e}")))?;
        for line in &bed_lines[start..end] {
            writeln!(bed_file, "{line}")
                .map_err(|e| PyIOError::new_err(format!("failed to write chunk BED: {e}")))?;
        }

        chunks.push(ChunkTask {
            index: i,
            bed_path: bed_chunk_path,
            stdout_path,
            stderr_path,
        });
        start = end;
    }

    let timeout = Duration::from_secs(3600);
    let next_index = Arc::new(Mutex::new(0usize));
    let (sender, receiver) = channel::<ChunkOutcome>();
    let chunks_ref = &chunks;

    thread::scope(|scope| {
        let num_workers = n_threads.min(chunks_ref.len());
        for _ in 0..num_workers {
            let sender: Sender<ChunkOutcome> = sender.clone();
            let next_index = Arc::clone(&next_index);
            scope.spawn(move || {
                loop {
                    let idx = {
                        let mut guard = next_index.lock().unwrap();
                        let idx = *guard;
                        if idx >= chunks_ref.len() {
                            break;
                        }
                        *guard += 1;
                        idx
                    };
                    let task = &chunks_ref[idx];
                    let outcome = run_freebayes_chunk(
                        freebayes_path,
                        cram_path,
                        fasta_path,
                        options,
                        task,
                        timeout,
                    );
                    if sender.send(outcome).is_err() {
                        break;
                    }
                }
            });
        }
    });

    let mut outcomes: HashMap<usize, ChunkOutcome> = HashMap::with_capacity(num_chunks);
    for _ in 0..num_chunks {
        let outcome = receiver
            .recv()
            .map_err(|e| FreebayesError::new_err(format!("worker channel error: {e}")))?;
        outcomes.insert(outcome.index(), outcome);
    }

    let mut stdout_paths: Vec<std::path::PathBuf> = Vec::with_capacity(num_chunks);
    for i in 0..num_chunks {
        match outcomes.remove(&i).expect("missing chunk outcome") {
            ChunkOutcome::Success { stdout_path, .. } => stdout_paths.push(stdout_path),
            ChunkOutcome::Failure {
                index,
                exit_code,
                stderr,
            } => {
                return Err(FreebayesError::new_err(format!(
                    "freebayes failed for chunk {index} (exit code {exit_code:?}): {stderr}"
                )));
            }
        }
    }

    let merged = merge_vcf_chunks(&stdout_paths)?;

    if let Some(path) = output_path {
        let mut out = File::create(path)
            .map_err(|e| PyIOError::new_err(format!("failed to create output VCF {path}: {e}")))?;
        out.write_all(merged.as_bytes())
            .map_err(|e| PyIOError::new_err(format!("failed to write output VCF {path}: {e}")))?;
    }

    Ok(merged)
}

impl ChunkOutcome {
    fn index(&self) -> usize {
        match self {
            ChunkOutcome::Success { index, .. } => *index,
            ChunkOutcome::Failure { index, .. } => *index,
        }
    }
}
