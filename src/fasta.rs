use anyhow::{Context, Result};
use fastx::FastX::{fasta_iter, reader_from_path, FastXRead};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Cursor};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Chromosome
{
    pub name: String,
    pub length: usize,
    pub sequence: Vec<u8>,
}

#[derive(Clone)]
pub struct Genome
{
    pub chromosomes: HashMap<String, Chromosome>,
}

impl Genome
{
    pub fn from_bytes(data: Vec<u8>, is_gzipped: bool) -> Result<Self>
    {
        let reader: Box<dyn BufRead> = if is_gzipped
        {
            Box::new(BufReader::new(MultiGzDecoder::new(Cursor::new(data))))
        }
        else
        {
            Box::new(BufReader::new(Cursor::new(data)))
        };

        Self::parse_fasta(reader)
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(path: &str) -> Result<Self>
    {
        if path.starts_with("http://") || path.starts_with("https://")
        {
            return Err(anyhow::anyhow!(
                "HTTP URLs not supported in native build. Please download the file first."
            ));
        }

        let reader = reader_from_path(Path::new(path))
            .with_context(|| format!("Failed to open FASTA file: {}", path))?;

        Self::parse_fasta(reader)
    }

    fn parse_fasta(reader: Box<dyn BufRead>) -> Result<Self>
    {
        let mut chromosomes = HashMap::new();

        for result in fasta_iter(reader)
        {
            let record = result.with_context(|| "Failed to parse FASTA record")?;

            let name = record.id().to_string();
            // Get sequence and convert to uppercase for consistency
            let sequence: Vec<u8> = record.seq().iter().map(|&b| b.to_ascii_uppercase()).collect();
            let length = sequence.len();

            let chr = Chromosome {
                name: name.clone(),
                length,
                sequence,
            };

            chromosomes.insert(name, chr);
        }

        if chromosomes.is_empty()
        {
            return Err(anyhow::anyhow!("No sequences found in FASTA file"));
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
