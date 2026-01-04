use anyhow::{Context, Result};
use fastx::indexed::IndexedFastXReader;
use fastx::FastX::{fasta_iter, reader_from_path, FastXRead};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor};
use std::path::Path;

#[cfg(not(target_arch = "wasm32"))]
use fastx::remote::RemoteReader;

/// Default block size for remote file caching (64KB as per fastx defaults)
#[cfg(not(target_arch = "wasm32"))]
const DEFAULT_REMOTE_BLOCK_SIZE: u64 = 64 * 1024;

/// Chromosome metadata (name and length)
#[derive(Debug, Clone)]
pub struct ChromosomeInfo
{
    pub name: String,
    pub length: usize,
}

/// Chromosome with full sequence data (for non-indexed mode)
#[derive(Debug, Clone)]
pub struct Chromosome
{
    pub name: String,
    pub length: usize,
    pub sequence: Vec<u8>,
}

/// Genome storage - either full sequences or indexed for lazy loading
pub enum GenomeData
{
    /// Full sequences loaded in memory (small files or non-indexed)
    Full(HashMap<String, Chromosome>),
    /// Indexed reader for lazy loading (bgzip with .fai + .gzi)
    IndexedLocal(IndexedFastXReader<File>),
    /// Remote indexed reader (HTTP/HTTPS with range requests) - native only
    #[cfg(not(target_arch = "wasm32"))]
    IndexedRemote(IndexedFastXReader<RemoteReader>),
}

/// Genome wrapper with chromosome metadata and optional lazy loading
pub struct Genome
{
    /// Chromosome metadata (always available)
    pub chromosomes: HashMap<String, ChromosomeInfo>,
    /// Genome data storage
    data: GenomeData,
    /// Cache for recently accessed chromosome data (for indexed mode)
    /// Stores full Chromosome structs for renderer compatibility
    chromosome_cache: HashMap<String, Chromosome>,
}

impl Clone for Genome
{
    fn clone(&self) -> Self
    {
        // For indexed readers, we only clone the metadata and cache
        // The caller should re-open indexed files if needed
        match &self.data
        {
            GenomeData::Full(chromosomes) =>
            {
                let full_chromosomes = chromosomes.clone();
                let chromosome_info = full_chromosomes
                    .iter()
                    .map(|(name, chr)| {
                        (
                            name.clone(),
                            ChromosomeInfo {
                                name: chr.name.clone(),
                                length: chr.length,
                            },
                        )
                    })
                    .collect();
                Self {
                    chromosomes: chromosome_info,
                    data: GenomeData::Full(full_chromosomes),
                    chromosome_cache: self.chromosome_cache.clone(),
                }
            }
            _ =>
            {
                // For indexed modes, we can't clone the reader
                // Return a genome with just metadata - caller should re-open if needed
                panic!("Cannot clone indexed genome - re-open the file instead");
            }
        }
    }
}

impl Genome
{
    /// Load genome from bytes (for WASM or in-memory data)
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

        Self::parse_fasta_full(reader)
    }

    /// Load genome from a local file path
    /// Automatically uses indexed loading if .fai and .gzi files are available
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(path: &str) -> Result<Self>
    {
        if path.starts_with("http://") || path.starts_with("https://")
        {
            return Self::from_url(path);
        }

        let path_obj = Path::new(path);

        // Try indexed loading first (for bgzip-compressed files with indexes)
        if let Ok(genome) = Self::try_indexed_local(path_obj)
        {
            return Ok(genome);
        }

        // Fall back to full loading
        let reader = reader_from_path(path_obj)
            .with_context(|| format!("Failed to open FASTA file: {}", path))?;

        Self::parse_fasta_full(reader)
    }

    /// Try to open as an indexed local file
    #[cfg(not(target_arch = "wasm32"))]
    fn try_indexed_local(path: &Path) -> Result<Self>
    {
        let reader = IndexedFastXReader::from_path(path)
            .with_context(|| format!("Failed to open indexed FASTA: {}", path.display()))?;

        let chromosomes = reader
            .index()
            .entries()
            .map(|entry| {
                (
                    entry.name.clone(),
                    ChromosomeInfo {
                        name: entry.name.clone(),
                        length: entry.length as usize,
                    },
                )
            })
            .collect();

        Ok(Self {
            chromosomes,
            data: GenomeData::IndexedLocal(reader),
            chromosome_cache: HashMap::new(),
        })
    }

    /// Load genome from a remote URL
    /// Requires the URL to point to a bgzip-compressed file with .fai and .gzi indexes
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_url(url: &str) -> Result<Self>
    {
        // Derive index URLs from the data URL
        let fai_url = format!("{}.fai", url);
        let gzi_url = format!("{}.gzi", url);

        Self::from_url_with_indexes(url, &fai_url, &gzi_url)
    }

    /// Load genome from a remote URL with explicit index URLs
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_url_with_indexes(data_url: &str, fai_url: &str, gzi_url: &str) -> Result<Self>
    {
        let reader = IndexedFastXReader::from_url(data_url, fai_url, gzi_url)
            .with_context(|| format!("Failed to open remote indexed FASTA: {}", data_url))?;

        let chromosomes = reader
            .index()
            .entries()
            .map(|entry| {
                (
                    entry.name.clone(),
                    ChromosomeInfo {
                        name: entry.name.clone(),
                        length: entry.length as usize,
                    },
                )
            })
            .collect();

        Ok(Self {
            chromosomes,
            data: GenomeData::IndexedRemote(reader),
            chromosome_cache: HashMap::new(),
        })
    }

    /// Load genome from a remote URL without requiring indexes
    /// This fetches the entire file into memory - suitable for smaller remote files
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_url_unindexed(url: &str, block_size: Option<u64>) -> Result<Self>
    {
        // Create a remote reader with optional custom block size
        let remote_reader = match block_size
        {
            Some(size) => RemoteReader::new(url)
                .map(|r| r.with_block_size(size))
                .with_context(|| format!("Failed to create remote reader for: {}", url))?,
            None => RemoteReader::new(url)
                .with_context(|| format!("Failed to create remote reader for: {}", url))?,
        };

        // Wrap in BufReader and use the standard FASTA iterator
        let buffered_reader = BufReader::new(remote_reader);
        Self::parse_fasta_full(Box::new(buffered_reader))
    }

    /// Parse FASTA and load all sequences into memory
    fn parse_fasta_full(reader: Box<dyn BufRead>) -> Result<Self>
    {
        let mut full_chromosomes = HashMap::new();

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

            full_chromosomes.insert(name, chr);
        }

        if full_chromosomes.is_empty()
        {
            return Err(anyhow::anyhow!("No sequences found in FASTA file"));
        }

        let chromosome_info = full_chromosomes
            .iter()
            .map(|(name, chr)| {
                (
                    name.clone(),
                    ChromosomeInfo {
                        name: chr.name.clone(),
                        length: chr.length,
                    },
                )
            })
            .collect();

        Ok(Genome {
            chromosomes: chromosome_info,
            data: GenomeData::Full(full_chromosomes),
            chromosome_cache: HashMap::new(),
        })
    }

    /// Get a chromosome by name (returns reference to ChromosomeInfo)
    pub fn get_chromosome_info(&self, name: &str) -> Option<&ChromosomeInfo>
    {
        self.chromosomes.get(name)
    }

    /// Check if this genome uses indexed (lazy) loading
    pub fn is_indexed(&self) -> bool
    {
        !matches!(self.data, GenomeData::Full(_))
    }

    /// Get a chromosome by name for rendering
    /// For full-loaded genomes, returns a reference directly
    /// For indexed genomes, fetches and caches the chromosome first
    pub fn get_chromosome(&self, chr_name: &str) -> Option<&Chromosome>
    {
        match &self.data
        {
            GenomeData::Full(chromosomes) => chromosomes.get(chr_name),
            _ =>
            {
                // For indexed genomes, check cache
                self.chromosome_cache.get(chr_name)
            }
        }
    }

    /// Ensure a chromosome is loaded and cached (for indexed genomes)
    /// Call this before get_chromosome() when using indexed mode
    pub fn ensure_chromosome_loaded(&mut self, chr_name: &str) -> Result<()>
    {
        match &mut self.data
        {
            GenomeData::Full(_) =>
            {
                // Already loaded
                Ok(())
            }
            GenomeData::IndexedLocal(reader) =>
            {
                if !self.chromosome_cache.contains_key(chr_name)
                {
                    let info = self
                        .chromosomes
                        .get(chr_name)
                        .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not found", chr_name))?;
                    let length = info.length;

                    let record = reader
                        .fetch(chr_name)
                        .with_context(|| format!("Failed to fetch chromosome '{}'", chr_name))?;
                    let sequence: Vec<u8> =
                        record.seq().iter().map(|&b| b.to_ascii_uppercase()).collect();

                    let chr = Chromosome {
                        name: chr_name.to_string(),
                        length,
                        sequence,
                    };
                    self.chromosome_cache.insert(chr_name.to_string(), chr);
                }
                Ok(())
            }
            #[cfg(not(target_arch = "wasm32"))]
            GenomeData::IndexedRemote(reader) =>
            {
                if !self.chromosome_cache.contains_key(chr_name)
                {
                    let info = self
                        .chromosomes
                        .get(chr_name)
                        .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not found", chr_name))?;
                    let length = info.length;

                    let record = reader
                        .fetch(chr_name)
                        .with_context(|| format!("Failed to fetch chromosome '{}'", chr_name))?;
                    let sequence: Vec<u8> =
                        record.seq().iter().map(|&b| b.to_ascii_uppercase()).collect();

                    let chr = Chromosome {
                        name: chr_name.to_string(),
                        length,
                        sequence,
                    };
                    self.chromosome_cache.insert(chr_name.to_string(), chr);
                }
                Ok(())
            }
        }
    }

    /// Get the full sequence for a chromosome
    /// For indexed genomes, this fetches and caches the sequence
    pub fn get_full_sequence(&mut self, chr_name: &str) -> Result<&[u8]>
    {
        self.ensure_chromosome_loaded(chr_name)?;
        match &self.data
        {
            GenomeData::Full(chromosomes) =>
            {
                chromosomes
                    .get(chr_name)
                    .map(|c| c.sequence.as_slice())
                    .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not found", chr_name))
            }
            _ => self
                .chromosome_cache
                .get(chr_name)
                .map(|c| c.sequence.as_slice())
                .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not cached", chr_name)),
        }
    }

    /// Get a sequence range for a chromosome
    /// For indexed genomes, uses efficient range fetching
    pub fn get_sequence_range(
        &mut self,
        chr_name: &str,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>>
    {
        match &mut self.data
        {
            GenomeData::Full(chromosomes) =>
            {
                let chr = chromosomes
                    .get(chr_name)
                    .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not found", chr_name))?;
                let start = start.min(chr.length);
                let end = end.min(chr.length);
                Ok(chr.sequence[start..end].to_vec())
            }
            GenomeData::IndexedLocal(reader) =>
            {
                // Check cache first
                if let Some(cached) = self.chromosome_cache.get(chr_name)
                {
                    let start = start.min(cached.length);
                    let end = end.min(cached.length);
                    return Ok(cached.sequence[start..end].to_vec());
                }

                // Fetch the range directly
                let data = reader
                    .fetch_range(chr_name, start as u64, end as u64)
                    .with_context(|| {
                        format!("Failed to fetch range {}:{}-{}", chr_name, start, end)
                    })?;
                // Convert to uppercase
                Ok(data.iter().map(|&b| b.to_ascii_uppercase()).collect())
            }
            #[cfg(not(target_arch = "wasm32"))]
            GenomeData::IndexedRemote(reader) =>
            {
                // Check cache first
                if let Some(cached) = self.chromosome_cache.get(chr_name)
                {
                    let start = start.min(cached.length);
                    let end = end.min(cached.length);
                    return Ok(cached.sequence[start..end].to_vec());
                }

                // Fetch the range directly
                let data = reader
                    .fetch_range(chr_name, start as u64, end as u64)
                    .with_context(|| {
                        format!("Failed to fetch range {}:{}-{}", chr_name, start, end)
                    })?;
                // Convert to uppercase
                Ok(data.iter().map(|&b| b.to_ascii_uppercase()).collect())
            }
        }
    }

    /// Clear the chromosome cache (to free memory)
    pub fn clear_cache(&mut self)
    {
        self.chromosome_cache.clear();
    }

    /// Get sequence region for a chromosome (compatibility method)
    pub fn get_sequence(&self, chr: &str, start: usize, end: usize) -> Option<&[u8]>
    {
        match &self.data
        {
            GenomeData::Full(chromosomes) => chromosomes
                .get(chr)
                .and_then(|c| c.sequence.get(start..end)),
            _ =>
            {
                // For indexed genomes, check cache only
                self.chromosome_cache
                    .get(chr)
                    .and_then(|c| c.sequence.get(start..end))
            }
        }
    }

    /// Check if a sequence exists in the genome
    /// For indexed genomes, this checks the index without fetching the sequence
    pub fn has_sequence(&self, seq_id: &str) -> bool
    {
        self.chromosomes.contains_key(seq_id)
    }

    /// Get all chromosome/sequence names in the genome
    pub fn sequence_names(&self) -> Vec<String>
    {
        self.chromosomes.keys().cloned().collect()
    }

    /// Get the total number of sequences in the genome
    pub fn sequence_count(&self) -> usize
    {
        self.chromosomes.len()
    }
}

/// Helper struct for renderer compatibility - wraps either a full Chromosome or cached data
pub struct ChromosomeView<'a>
{
    pub name: &'a str,
    pub length: usize,
    sequence: ChromosomeSequence<'a>,
}

enum ChromosomeSequence<'a>
{
    Full(&'a [u8]),
    Cached(&'a [u8]),
}

impl<'a> ChromosomeView<'a>
{
    pub fn get_gc_content(&self, start: usize, end: usize) -> f32
    {
        let seq = match &self.sequence
        {
            ChromosomeSequence::Full(s) => *s,
            ChromosomeSequence::Cached(s) => *s,
        };

        let start = start.min(seq.len());
        let end = end.min(seq.len());

        if start >= end
        {
            return 0.0;
        }

        let region = &seq[start..end];
        let gc_count = region.iter().filter(|&&b| b == b'G' || b == b'C').count();

        gc_count as f32 / region.len() as f32
    }

    pub fn get_base_at(&self, pos: usize) -> Option<char>
    {
        let seq = match &self.sequence
        {
            ChromosomeSequence::Full(s) => *s,
            ChromosomeSequence::Cached(s) => *s,
        };
        seq.get(pos).map(|&b| b as char)
    }

    pub fn sequence(&self) -> &[u8]
    {
        match &self.sequence
        {
            ChromosomeSequence::Full(s) => s,
            ChromosomeSequence::Cached(s) => s,
        }
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
