use anyhow::{Context, Result};
use memmap2::Mmap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct Chromosome
{
    pub name: String,
    pub length: usize,
    pub sequence: Vec<u8>,
}

pub struct Genome
{
    pub chromosomes: HashMap<String, Chromosome>,
}

impl Genome
{
    pub fn from_file(path: &str) -> Result<Self>
    {
        let file =
            File::open(path).with_context(|| format!("Failed to open FASTA file: {}", path))?;

        let reader = BufReader::new(file);
        let mut chromosomes = HashMap::new();
        let mut current_name = String::new();
        let mut current_seq = Vec::new();

        for line in reader.lines()
        {
            let line = line?;
            let line = line.trim();

            if line.is_empty()
            {
                continue;
            }

            if line.starts_with('>')
            {
                // Save previous chromosome
                if !current_name.is_empty()
                {
                    let chr = Chromosome {
                        name: current_name.clone(),
                        length: current_seq.len(),
                        sequence: current_seq.clone(),
                    };
                    chromosomes.insert(current_name.clone(), chr);
                }

                // Start new chromosome
                current_name = line[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                current_seq = Vec::new();
            }
            else
            {
                // Add sequence data (convert to uppercase for consistency)
                current_seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
            }
        }

        // Save last chromosome
        if !current_name.is_empty()
        {
            let chr = Chromosome {
                name: current_name.clone(),
                length: current_seq.len(),
                sequence: current_seq,
            };
            chromosomes.insert(current_name, chr);
        }

        Ok(Genome { chromosomes })
    }

    pub fn get_sequence(&self, chr: &str, start: usize, end: usize) -> Option<&[u8]>
    {
        self.chromosomes
            .get(chr)
            .and_then(|c| c.sequence.get(start..end))
    }
}

impl Chromosome
{
    pub fn get_gc_content(&self, start: usize, end: usize) -> f32
    {
        let start = start.min(self.length);
        let end = end.min(self.length);

        if start >= end
        {
            return 0.0;
        }

        let seq = &self.sequence[start..end];
        let gc_count = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();

        gc_count as f32 / seq.len() as f32
    }

    pub fn get_base_at(&self, pos: usize) -> Option<char>
    {
        self.sequence.get(pos).map(|&b| b as char)
    }
}
