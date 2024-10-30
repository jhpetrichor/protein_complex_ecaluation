use std::path::{Path, PathBuf};

#[derive(Debug)]
pub struct Options {
    pub ppi_path: PathBuf,
    pub ref_complex_path: PathBuf,
    pub pre_complex_path: PathBuf,
    pub min_size: u32,
    pub threshold: f64,
}
