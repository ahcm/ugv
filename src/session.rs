use anyhow::Result;
use serde::{Deserialize, Serialize};

/// Session data to persist viewer state
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Session
{
    pub fasta_path: String,
    pub gff_path: String,
    pub bam_path: String,
    pub selected_chromosome: Option<String>,
    pub viewport_start: usize,
    pub viewport_end: usize,
    pub show_chromosome_panel: bool,
    pub show_amino_acids: bool,
    pub show_coverage: bool,
    pub show_alignments: bool,
    pub show_variants: bool,
    pub chromosome_sort: String, // "Natural", "Alphabetical", or "Size"
}

impl Session
{
    /// Create a new session with default values
    pub fn new() -> Self
    {
        Self {
            fasta_path: String::new(),
            gff_path: String::new(),
            bam_path: String::new(),
            selected_chromosome: None,
            viewport_start: 0,
            viewport_end: 1000000,
            show_chromosome_panel: true,
            show_amino_acids: false,
            show_coverage: true,
            show_alignments: true,
            show_variants: false,
            chromosome_sort: "Natural".to_string(),
        }
    }

    /// Save session to JSON string
    pub fn to_json(&self) -> Result<String>
    {
        Ok(serde_json::to_string_pretty(self)?)
    }

    /// Load session from JSON string
    pub fn from_json(json: &str) -> Result<Self>
    {
        Ok(serde_json::from_str(json)?)
    }

    /// Save session to file (native only)
    #[cfg(not(target_arch = "wasm32"))]
    pub fn save_to_file(&self, path: &str) -> Result<()>
    {
        use std::fs;
        let json = self.to_json()?;
        fs::write(path, json)?;
        Ok(())
    }

    /// Load session from file (native only)
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load_from_file(path: &str) -> Result<Self>
    {
        use std::fs;
        let json = fs::read_to_string(path)?;
        Self::from_json(&json)
    }
}

impl Default for Session
{
    fn default() -> Self
    {
        Self::new()
    }
}
