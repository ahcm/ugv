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

/// FAI index entry for a sequence
#[derive(Debug, Clone)]
pub struct FaiEntry
{
    pub name: String,
    pub length: u64,
    pub offset: u64,
    pub line_bases: u64,
    pub line_bytes: u64,
}

/// GZip index entry for a BGZF block
#[derive(Debug, Clone)]
pub struct GziEntry
{
    pub compressed_offset: u64,
    pub uncompressed_offset: u64,
    pub compressed_size: u64,
    pub uncompressed_size: u64,
}

/// WASM remote indexed FASTA data (stores URLs and index data for lazy loading)
#[cfg(target_arch = "wasm32")]
#[derive(Clone)]
pub struct WasmRemoteIndexed
{
    pub data_url: String,
    pub fai_entries: Vec<FaiEntry>,
    pub gzi_entries: Vec<GziEntry>,
}

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
    /// WASM remote indexed FASTA (stores URLs and index data for lazy loading via fetch)
    #[cfg(target_arch = "wasm32")]
    WasmRemoteIndexed(WasmRemoteIndexed),
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

    /// Load genome from a URL for WASM, checking for index files
    /// This function is async-compatible and returns a genome that can fetch ranges
    /// Supports both uncompressed FASTA and BGZF-compressed FASTA with .fai and .gzi indexes
    #[cfg(target_arch = "wasm32")]
    pub fn from_url_wasm(url: &str, fai_data: Vec<u8>, gzi_data: Option<Vec<u8>>) -> Result<Self>
    {
        // Parse FAI index
        let fai_entries = Self::parse_fai(&fai_data)?;

        // Parse GZI index if provided (for gzipped files)
        let gzi_entries = if let Some(gzi) = gzi_data
        {
            if url.ends_with(".gz")
            {
                Some(Self::parse_gzi(&gzi)?)
            }
            else
            {
                None
            }
        }
        else
        {
            None
        };

        // For gzipped files without GZI, we can't do efficient range fetching
        if url.ends_with(".gz") && gzi_entries.is_none()
        {
            return Err(anyhow::anyhow!(
                "Gzipped file requires .gzi index for efficient range fetching"
            ));
        }

        let chromosomes = fai_entries
            .iter()
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

        let wasm_indexed = WasmRemoteIndexed {
            data_url: url.to_string(),
            fai_entries,
            gzi_entries: gzi_entries.unwrap_or_default(),
        };

        Ok(Self {
            chromosomes,
            data: GenomeData::WasmRemoteIndexed(wasm_indexed),
            chromosome_cache: HashMap::new(),
        })
    }

    /// Parse FAI (FASTA index) file format
    /// Format: sequence_name length offset line_bases line_bytes
    #[cfg(target_arch = "wasm32")]
    fn parse_fai(data: &[u8]) -> Result<Vec<FaiEntry>>
    {
        let text = std::str::from_utf8(data)
            .with_context(|| "FAI data is not valid UTF-8")?;

        let mut entries = Vec::new();
        for line in text.lines()
        {
            let line = line.trim();
            if line.is_empty()
            {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() != 5
            {
                return Err(anyhow::anyhow!("Invalid FAI line: {}", line));
            }

            entries.push(FaiEntry {
                name: parts[0].to_string(),
                length: parts[1].parse()
                    .with_context(|| format!("Invalid length in FAI: {}", parts[1]))?,
                offset: parts[2].parse()
                    .with_context(|| format!("Invalid offset in FAI: {}", parts[2]))?,
                line_bases: parts[3].parse()
                    .with_context(|| format!("Invalid line_bases in FAI: {}", parts[3]))?,
                line_bytes: parts[4].parse()
                    .with_context(|| format!("Invalid line_bytes in FAI: {}", parts[4]))?,
            });
        }
        Ok(entries)
    }

    /// Parse GZI (GZIP index) file format
    /// Binary format: 8 bytes (number of offsets) + N * 16 bytes (offset entries)
    /// Each entry: 8 bytes compressed offset + 8 bytes uncompressed offset
    #[cfg(target_arch = "wasm32")]
    fn parse_gzi(data: &[u8]) -> Result<Vec<GziEntry>>
    {
        if data.len() < 8
        {
            return Err(anyhow::anyhow!("Invalid GZI file: too small"));
        }

        // First 8 bytes: number of offsets (little-endian u64)
        let num_offsets = u64::from_le_bytes(
            data[0..8].try_into().unwrap()
        ) as usize;

        let expected_size = 8 + (num_offsets * 16);
        if data.len() < expected_size
        {
            return Err(anyhow::anyhow!(
                "Invalid GZI file: expected {} bytes, got {}", expected_size, data.len()
            ));
        }

        let mut entries = Vec::new();
        for i in 0..num_offsets
        {
            let offset = 8 + (i * 16);
            let compressed_offset = u64::from_le_bytes(
                data[offset..offset + 8].try_into().unwrap()
            );
            let uncompressed_offset = u64::from_le_bytes(
                data[offset + 8..offset + 16].try_into().unwrap()
            );

            // Calculate sizes (difference from next entry, or 0 for last)
            let (compressed_size, uncompressed_size) = if i + 1 < num_offsets
            {
                let next_compressed = u64::from_le_bytes(
                    data[offset + 16..offset + 24].try_into().unwrap()
                );
                let next_uncompressed = u64::from_le_bytes(
                    data[offset + 24..offset + 32].try_into().unwrap()
                );
                (next_compressed - compressed_offset, next_uncompressed - uncompressed_offset)
            }
            else
            {
                (0, 0)
            };

            entries.push(GziEntry {
                compressed_offset,
                uncompressed_offset,
                compressed_size,
                uncompressed_size,
            });
        }
        Ok(entries)
    }

    /// Fetch a byte range from a URL using HTTP range request
    /// This is used internally by WASM indexed genome
    #[cfg(target_arch = "wasm32")]
    pub async fn fetch_url_range(url: &str, start: u64, end: u64) -> Result<Vec<u8>>
    {
        use wasm_bindgen::JsCast;
        use wasm_bindgen_futures::JsFuture;
        use web_sys::{RequestInit, RequestMode, Response};

        let window = web_sys::window().ok_or_else(|| anyhow::anyhow!("No window object"))?;

        // Build URL with range parameter (since we can't easily set Range header via fetch API)
        // We'll use the full fetch and let the server handle the range
        let range_header = format!("bytes={}-{}", start, end.saturating_sub(1));

        let mut init = RequestInit::new();
        init.set_method("GET");
        init.set_mode(RequestMode::Cors);

        // Set Range header using the headers() method
        let headers_val = js_sys::Object::new();
        js_sys::Reflect::set(&headers_val, &"Range".into(), &range_header.into()).unwrap();
        init.set_headers(&headers_val.into());

        let response_promise: js_sys::Promise = window.fetch_with_str_and_init(url, &init).into();
        let response_value = JsFuture::from(response_promise)
            .await
            .map_err(|e| anyhow::anyhow!("Fetch failed: {:?}", e))?;

        let response: Response = response_value
            .dyn_into()
            .map_err(|_| anyhow::anyhow!("Response cast failed"))?;

        // Accept both 200 and 206 (Partial Content)
        let status = response.status();
        if status != 200 && status != 206
        {
            return Err(anyhow::anyhow!(
                "HTTP error {}: {}", status, response.status_text()
            ));
        }

        let array_buffer_promise = response
            .array_buffer()
            .map_err(|_| anyhow::anyhow!("Failed to get array buffer"))?;

        let array_buffer_value = JsFuture::from(array_buffer_promise)
            .await
            .map_err(|e| anyhow::anyhow!("Failed to read array buffer: {:?}", e))?;

        let uint8_array = js_sys::Uint8Array::new(&array_buffer_value);
        Ok(uint8_array.to_vec())
    }

    /// Async version of chromosome loading for WASM indexed genomes
    /// Fetches a chromosome from the remote URL using the index data
    /// Handles both uncompressed and BGZF-compressed FASTA files
    #[cfg(target_arch = "wasm32")]
    pub async fn load_chromosome_async(&mut self, chr_name: &str) -> Result<()>
    {
        // If already cached, return immediately
        if self.chromosome_cache.contains_key(chr_name)
        {
            return Ok(());
        }

        let (data_url, fai_entries, gzi_entries, is_gzipped) = match &self.data
        {
            GenomeData::WasmRemoteIndexed(wasm_indexed) => {
                let is_gzipped = wasm_indexed.data_url.ends_with(".gz");
                (
                    wasm_indexed.data_url.clone(),
                    wasm_indexed.fai_entries.clone(),
                    wasm_indexed.gzi_entries.clone(),
                    is_gzipped,
                )
            }
            _ => return Ok(()), // Not a WASM indexed genome
        };

        let info = self
            .chromosomes
            .get(chr_name)
            .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not found", chr_name))?;

        let length = info.length;

        // Find the FAI entry for this chromosome
        let fai_entry = fai_entries
            .iter()
            .find(|e| e.name == chr_name)
            .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not in FAI index", chr_name))?;

        let sequence = if is_gzipped
        {
            // For BGZF-compressed files, use the GZI index to fetch and decompress blocks
            Self::fetch_bgzf_sequence(&data_url, fai_entry, &gzi_entries, length).await?
        }
        else
        {
            // For uncompressed files, fetch the range directly
            let data_end = fai_entry.offset + (length as u64);
            let raw_data = Self::fetch_url_range(&data_url, fai_entry.offset, data_end).await?;
            Self::parse_fasta_sequence(&raw_data)
        };

        if sequence.len() != length
        {
            return Err(anyhow::anyhow!(
                "Sequence length mismatch for '{}': expected {}, got {}",
                chr_name, length, sequence.len()
            ));
        }

        let chr = Chromosome {
            name: chr_name.to_string(),
            length,
            sequence,
        };
        self.chromosome_cache.insert(chr_name.to_string(), chr);
        Ok(())
    }

    /// Fetch and decompress a sequence from a BGZF-compressed file using GZI index
    #[cfg(target_arch = "wasm32")]
    async fn fetch_bgzf_sequence(
        url: &str,
        fai_entry: &FaiEntry,
        gzi_entries: &[GziEntry],
        expected_length: usize,
    ) -> Result<Vec<u8>>
    {
        // Calculate the byte range we need in uncompressed space
        let start_uncompressed = fai_entry.offset;
        let end_uncompressed = fai_entry.offset + fai_entry.length;

        // Find which GZI blocks overlap with our range
        let start_block_idx = gzi_entries
            .partition_point(|e| e.uncompressed_offset < start_uncompressed)
            .saturating_sub(1);

        let end_block_idx = gzi_entries
            .partition_point(|e| e.uncompressed_offset < end_uncompressed)
            .min(gzi_entries.len());

        if start_block_idx >= gzi_entries.len()
        {
            return Err(anyhow::anyhow!("No GZI blocks found for sequence"));
        }

        // Fetch all required blocks (in one combined range request if possible)
        let first_block = &gzi_entries[start_block_idx];
        let last_block = &gzi_entries[end_block_idx.saturating_sub(1)];

        // Fetch from first block's compressed offset to end of last block
        let fetch_start = first_block.compressed_offset;
        // Estimate the end: last block offset + 64KB (max BGZF block size)
        let fetch_end = last_block.compressed_offset + (64 * 1024);

        let compressed_data = Self::fetch_url_range(url, fetch_start, fetch_end).await?;

        // Decompress blocks and extract our range
        let mut result = Vec::with_capacity(expected_length);
        let mut current_uncompressed_offset = 0u64;

        for block_idx in start_block_idx..end_block_idx
        {
            let block = &gzi_entries[block_idx];

            // Calculate where this block is within our fetched data
            let block_start_in_fetch = (block.compressed_offset - fetch_start) as usize;
            let block_end_in_fetch = if block_idx + 1 < gzi_entries.len()
            {
                let next_block = &gzi_entries[block_idx + 1];
                (next_block.compressed_offset - fetch_start) as usize
            }
            else
            {
                compressed_data.len()
            };

            if block_end_in_fetch > compressed_data.len()
            {
                // Need to fetch more data
                return Err(anyhow::anyhow!("Block extends beyond fetched data"));
            }

            let block_data = &compressed_data[block_start_in_fetch..block_end_in_fetch];

            // Decompress this BGZF block
            let decompressed = Self::decompress_bgzf_block(block_data)?;

            // Calculate the offset within this decompressed block
            let block_start_offset = if block.uncompressed_offset >= start_uncompressed
            {
                0
            }
            else
            {
                (start_uncompressed - block.uncompressed_offset) as usize
            };

            let block_end_offset = if block.uncompressed_offset + decompressed.len() as u64 <= end_uncompressed
            {
                decompressed.len()
            }
            else
            {
                (end_uncompressed - block.uncompressed_offset) as usize
            };

            // Extract our range from this block
            if block_start_offset < decompressed.len() && block_end_offset > block_start_offset
            {
                let extract_end = block_end_offset.min(decompressed.len());
                result.extend_from_slice(&decompressed[block_start_offset..extract_end]);
            }

            current_uncompressed_offset = block.uncompressed_offset + decompressed.len() as u64;
        }

        // Parse the FASTA sequence from the decompressed data
        Ok(Self::parse_fasta_sequence(&result))
    }

    /// Decompress a single BGZF block
    /// BGZF blocks are valid gzip streams but with specific compression parameters
    #[cfg(target_arch = "wasm32")]
    fn decompress_bgzf_block(data: &[u8]) -> Result<Vec<u8>>
    {
        use flate2::read::GzDecoder;

        // Check for BGZF magic bytes (first 2 bytes should be 0x1f 0x8b)
        if data.len() < 18
        {
            return Err(anyhow::anyhow!("BGZF block too small: {}", data.len()));
        }

        // BGZF-specific: check XFL field at offset 8 (should be 4 for BGZF)
        // and extra flags at offset 3 (should be 4 for BGZF)
        if data[0] != 0x1f || data[1] != 0x8b
        {
            return Err(anyhow::anyhow!("Invalid gzip header"));
        }

        // Use flate2 to decompress
        let mut decoder = GzDecoder::new(data);
        let mut decompressed = Vec::new();
        std::io::Read::read_to_end(&mut decoder, &mut decompressed)
            .map_err(|e| anyhow::anyhow!("Decompression failed: {}", e))?;

        Ok(decompressed)
    }

    /// Parse FASTA sequence from raw bytes (removes newlines and header)
    #[cfg(target_arch = "wasm32")]
    fn parse_fasta_sequence(data: &[u8]) -> Vec<u8>
    {
        let mut result = Vec::new();
        let mut in_sequence = false;

        for line in data.split(|&b| b == b'\n')
        {
            if line.is_empty()
            {
                continue;
            }
            // Skip header line
            if line[0] == b'>'
            {
                in_sequence = true;
                continue;
            }
            // Add sequence bytes, converting to uppercase
            if in_sequence
            {
                for &byte in line
                {
                    if byte != b'\r' && byte != b'\n' && !byte.is_ascii_whitespace()
                    {
                        result.push(byte.to_ascii_uppercase());
                    }
                }
            }
        }
        result
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
            #[cfg(target_arch = "wasm32")]
            GenomeData::WasmRemoteIndexed(_) =>
            {
                // For WASM, chromosome loading is async and must be handled separately
                // Use ensure_chromosome_loaded_async instead
                Err(anyhow::anyhow!(
                    "WASM indexed genome requires async loading. Use load_chromosome_wasm_async() instead."
                ))
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
            #[cfg(target_arch = "wasm32")]
            GenomeData::WasmRemoteIndexed(_) =>
            {
                // For WASM indexed genomes, check cache first
                if let Some(cached) = self.chromosome_cache.get(chr_name)
                {
                    let start = start.min(cached.length);
                    let end = end.min(cached.length);
                    return Ok(cached.sequence[start..end].to_vec());
                }
                // Chromosome not loaded - needs async loading
                Err(anyhow::anyhow!(
                    "Chromosome '{}' not loaded. Call load_chromosome_async() first.",
                    chr_name
                ))
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
