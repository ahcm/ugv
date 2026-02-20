use anyhow::Result;
use serde::{Deserialize, Serialize};

fn default_tsv_label_font_size() -> f32
{
    11.0
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SessionTrackConfig
{
    pub track_type: String,
    pub height: f32,
    pub visible: bool,
    pub order: usize,
}

/// Session data to persist viewer state
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Session
{
    pub fasta_path: String,
    pub gff_path: String,
    pub bam_path: String,
    #[serde(default)]
    pub bam_path_2: String,
    #[serde(default)]
    pub tsv_path: String,
    pub selected_chromosome: Option<String>,
    pub viewport_start: usize,
    pub viewport_end: usize,
    pub show_chromosome_panel: bool,
    pub show_amino_acids: bool,
    pub show_coverage: bool,
    pub show_alignments: bool,
    pub show_variants: bool,
    #[serde(default)]
    pub track_order: Vec<String>,
    #[serde(default)]
    pub track_configs: Vec<SessionTrackConfig>,
    pub chromosome_sort: String, // "Natural", "Alphabetical", or "Size"
    #[serde(default = "default_max_reads_display")]
    pub max_reads_display: usize,
    #[serde(default = "default_bam_memory_budget_mib")]
    pub bam_memory_budget_mib: usize,
    #[serde(default = "default_tsv_label_font_size")]
    pub tsv_label_font_size: f32,
    #[serde(default)]
    pub tsv_label_font_monospace: bool,
}

fn default_max_reads_display() -> usize
{
    1000
}

fn default_bam_memory_budget_mib() -> usize
{
    4096
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
            bam_path_2: String::new(),
            tsv_path: String::new(),
            selected_chromosome: None,
            viewport_start: 0,
            viewport_end: 1000000,
            show_chromosome_panel: true,
            show_amino_acids: false,
            show_coverage: true,
            show_alignments: true,
            show_variants: false,
            track_order: Vec::new(),
            track_configs: Vec::new(),
            chromosome_sort: "Natural".to_string(),
            max_reads_display: 1000,
            bam_memory_budget_mib: default_bam_memory_budget_mib(),
            tsv_label_font_size: default_tsv_label_font_size(),
            tsv_label_font_monospace: false,
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
