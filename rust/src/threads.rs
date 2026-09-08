use std::env;

const DEFAULT_THREAD_CAP: usize = 8;
const ENV_VAR: &str = "METEOR_RUST_THREADS";

/// Resolve the number of threads to use for htslib CRAM decode.
///
/// - If `METEOR_RUST_THREADS` is unset, returns the smaller of 8 and the
///   number of logical CPUs reported by the system, defaulting to 1 when the
///   parallelism cannot be determined.
/// - If the variable is set, it is parsed as an unsigned integer and clamped
///   to the inclusive range `[1, available_parallelism()]`. Values that are
///   zero or exceed the number of logical CPUs fall back to a safe thread
///   count. A non-numeric value prints a warning to stderr before falling back.
pub(crate) fn resolve_thread_count() -> usize {
    match env::var(ENV_VAR) {
        Ok(value) => match value.parse::<usize>() {
            Ok(0) => default_thread_count(),
            Ok(n) => n.min(max_thread_count()),
            Err(_) => {
                let default = default_thread_count();
                eprintln!(
                    "METEOR_RUST_THREADS={value:?} is not an integer; using default thread count {default}"
                );
                default
            }
        },
        Err(_) => default_thread_count(),
    }
}

fn max_thread_count() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

fn default_thread_count() -> usize {
    DEFAULT_THREAD_CAP.min(max_thread_count())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_at_least_one() {
        // default_thread_count is independent of the environment and should
        // always return a positive number because available_parallelism >= 1.
        assert!(default_thread_count() >= 1);
    }
}
