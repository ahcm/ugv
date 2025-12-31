use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub enum Strand
{
    Forward,
    Reverse,
    Unknown,
}

#[derive(Debug, Clone)]
pub struct Feature
{
    pub seqid: String,
    pub source: String,
    pub feature_type: String,
    pub start: usize,
    pub end: usize,
    pub score: Option<f64>,
    pub strand: Strand,
    pub phase: Option<u8>,
    pub attributes: HashMap<String, String>,
}

impl Feature
{
    pub fn get_attribute(&self, key: &str) -> Option<&str>
    {
        self.attributes.get(key).map(|s| s.as_str())
    }

    pub fn name(&self) -> String
    {
        self.get_attribute("Name")
            .or_else(|| self.get_attribute("ID"))
            .or_else(|| self.get_attribute("gene"))
            .unwrap_or(&self.feature_type)
            .to_string()
    }

    pub fn color_by_type(&self) -> egui::Color32
    {
        match self.feature_type.as_str()
        {
            "gene" => egui::Color32::from_rgb(70, 130, 180),
            "mRNA" | "transcript" => egui::Color32::from_rgb(100, 150, 200),
            "exon" => egui::Color32::from_rgb(34, 139, 34),
            "CDS" => egui::Color32::from_rgb(255, 140, 0),
            "UTR" | "five_prime_UTR" | "three_prime_UTR" => egui::Color32::from_rgb(220, 220, 100),
            "intron" => egui::Color32::from_rgb(180, 180, 180),
            _ => egui::Color32::from_rgb(128, 128, 128),
        }
    }
}

pub fn parse_gff(path: &str) -> Result<Vec<Feature>>
{
    let file = File::open(path).with_context(|| format!("Failed to open GFF file: {}", path))?;

    let reader: Box<dyn BufRead> = if path.ends_with(".gz")
    {
        Box::new(BufReader::new(GzDecoder::new(file)))
    }
    else
    {
        Box::new(BufReader::new(file))
    };

    let mut features = Vec::new();

    for line in reader.lines()
    {
        let line = line?;
        let line = line.trim();

        // Skip comments and empty lines
        if line.is_empty() || line.starts_with('#')
        {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 8
        {
            continue;
        }

        let seqid = parts[0].to_string();
        let source = parts[1].to_string();
        let feature_type = parts[2].to_string();

        let start = parts[3]
            .parse::<usize>()
            .with_context(|| format!("Invalid start position: {}", parts[3]))?
            .saturating_sub(1); // Convert to 0-based

        let end = parts[4]
            .parse::<usize>()
            .with_context(|| format!("Invalid end position: {}", parts[4]))?;

        let score = if parts[5] == "."
        {
            None
        }
        else
        {
            Some(parts[5].parse::<f64>()?)
        };

        let strand = match parts[6]
        {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Unknown,
        };

        let phase = if parts[7] == "."
        {
            None
        }
        else
        {
            Some(parts[7].parse::<u8>()?)
        };

        let attributes = if parts.len() > 8
        {
            parse_attributes(parts[8])
        }
        else
        {
            HashMap::new()
        };

        features.push(Feature {
            seqid,
            source,
            feature_type,
            start,
            end,
            score,
            strand,
            phase,
            attributes,
        });
    }

    Ok(features)
}

fn parse_attributes(attr_str: &str) -> HashMap<String, String>
{
    let mut attributes = HashMap::new();

    for pair in attr_str.split(';')
    {
        let pair = pair.trim();
        if pair.is_empty()
        {
            continue;
        }

        if let Some((key, value)) = pair.split_once('=')
        {
            // GFF3 format
            attributes.insert(key.trim().to_string(), value.trim().trim_matches('"').to_string());
        }
        else if let Some((key, value)) = pair.split_once(' ')
        {
            // GTF format
            attributes.insert(key.trim().to_string(), value.trim().trim_matches('"').to_string());
        }
    }

    attributes
}
