use anyhow::{Context, Result};
use fastx::FastX::{FastXRead, fasta_iter};
#[cfg(target_arch = "wasm32")]
use fastx::fai::FaiEntry;

#[cfg(not(target_arch = "wasm32"))]
use fastx::FastX::reader_from_path;
#[cfg(not(target_arch = "wasm32"))]
use std::path::Path;

#[cfg(target_arch = "wasm32")]
#[derive(Clone, Debug)]
pub struct SimpleFaiIndex
{
    pub entries: HashMap<String, FaiEntry>,
}

#[cfg(target_arch = "wasm32")]
impl SimpleFaiIndex
{
    pub fn from_reader<R: BufRead>(mut reader: R) -> Result<Self>
    {
        let mut entries = HashMap::new();
        let mut line = String::new();
        while reader.read_line(&mut line)? > 0
        {
            let parts: Vec<&str> = line.trim().split('\t').collect();
            if parts.len() >= 5
            {
                let name = parts[0].to_string();
                let length = parts[1].parse()?;
                let offset = parts[2].parse()?;
                let line_bases = parts[3].parse()?;
                let line_width = parts[4].parse()?;

                // Construct FaiEntry manually since we can't depend on fastx constructors
                // Assuming fields are public, otherwise we might have issues.
                // If FaiEntry fields are private, we'll need to define our own FaiEntry too.
                let entry = FaiEntry {
                    name: name.clone(),
                    length,
                    offset,
                    line_bases,
                    line_width,
                };
                entries.insert(name, entry);
            }
            line.clear();
        }
        Ok(Self { entries })
    }

    pub fn get(&self, name: &str) -> Option<&FaiEntry>
    {
        self.entries.get(name)
    }
}

#[cfg(target_arch = "wasm32")]
#[derive(Clone, Debug)]
pub struct SimpleGziIndex
{
    entries: Vec<(u64, u64)>,
}

#[cfg(target_arch = "wasm32")]
impl SimpleGziIndex
{
    pub fn from_bytes(data: &[u8]) -> Result<Self>
    {
        use byteorder::{ByteOrder, LittleEndian};
        if data.len() < 8
        {
            return Err(anyhow::anyhow!("Invalid GZI data"));
        }
        let count = LittleEndian::read_u64(&data[0..8]);
        let mut entries = Vec::with_capacity(count as usize);
        let mut offset = 8;
        for _ in 0..count
        {
            if offset + 16 > data.len()
            {
                break;
            }
            let compressed = LittleEndian::read_u64(&data[offset..offset + 8]);
            let uncompressed = LittleEndian::read_u64(&data[offset + 8..offset + 16]);
            entries.push((compressed, uncompressed));
            offset += 16;
        }
        Ok(Self { entries })
    }

    pub fn entries(&self) -> &[(u64, u64)]
    {
        &self.entries
    }
}

use fastx::indexed::IndexedFastXReader;

use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor};

#[cfg(not(target_arch = "wasm32"))]
use fastx::remote::RemoteReader;

/// WASM remote indexed FASTA data (stores URLs and index data for lazy loading)
#[cfg(target_arch = "wasm32")]
#[derive(Clone)]
pub struct WasmRemoteIndexed
{
    pub data_url: String,
    pub fai_index: SimpleFaiIndex,
    pub gzi_index: Option<SimpleGziIndex>,
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
            #[cfg(target_arch = "wasm32")]
            GenomeData::WasmRemoteIndexed(wasm_indexed) => Self {
                chromosomes: self.chromosomes.clone(),
                data: GenomeData::WasmRemoteIndexed(wasm_indexed.clone()),
                chromosome_cache: self.chromosome_cache.clone(),
            },
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
        // Parse FAI index using fastx
        let fai_index = SimpleFaiIndex::from_reader(Cursor::new(fai_data))
            .map_err(|e| anyhow::anyhow!("Failed to parse FAI: {}", e))?;

        // Parse GZI index if provided (for gzipped files)
        let gzi_index = if let Some(gzi) = gzi_data
        {
            if url.ends_with(".gz")
            {
                Some(
                    SimpleGziIndex::from_bytes(&gzi)
                        .map_err(|e| anyhow::anyhow!("Failed to parse GZI: {}", e))?,
                )
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
        if url.ends_with(".gz") && gzi_index.is_none()
        {
            return Err(anyhow::anyhow!(
                "Gzipped file requires .gzi index for efficient range fetching"
            ));
        }

        let chromosomes = fai_index
            .entries
            .iter()
            .map(|(name, entry)| {
                (
                    name.clone(),
                    ChromosomeInfo {
                        name: name.clone(),
                        length: entry.length as usize,
                    },
                )
            })
            .collect();

        let wasm_indexed = WasmRemoteIndexed {
            data_url: url.to_string(),
            fai_index,
            gzi_index,
        };

        Ok(Self {
            chromosomes,
            data: GenomeData::WasmRemoteIndexed(wasm_indexed),
            chromosome_cache: HashMap::new(),
        })
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

        let range_header = format!("bytes={}-{}", start, end.saturating_sub(1));

        let init = RequestInit::new();
        init.set_method("GET");
        init.set_mode(RequestMode::Cors);

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

        let status = response.status();
        if status != 200 && status != 206
        {
            return Err(anyhow::anyhow!("HTTP error {}: {}", status, response.status_text()));
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

        let (data_url, fai_index, gzi_index, is_gzipped) = match &self.data
        {
            GenomeData::WasmRemoteIndexed(wasm_indexed) =>
            {
                let is_gzipped = wasm_indexed.data_url.ends_with(".gz");
                (
                    wasm_indexed.data_url.clone(),
                    wasm_indexed.fai_index.clone(),
                    wasm_indexed.gzi_index.clone(),
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
        let fai_entry = fai_index
            .get(chr_name)
            .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not in FAI index", chr_name))?;

        let sequence = if is_gzipped
        {
            // For BGZF-compressed files, use the GZI index to fetch and decompress blocks
            let gzi = gzi_index
                .ok_or_else(|| anyhow::anyhow!("Missing GZI index for compressed file"))?;
            Self::fetch_bgzf_sequence(&data_url, fai_entry, &gzi, length).await?
        }
        else
        {
            // For uncompressed files, fetch the range directly
            let data_end = fai_entry.offset_for_position(fai_entry.length as u64) as u64;
            let raw_data = Self::fetch_url_range(&data_url, fai_entry.offset, data_end).await?;
            Self::process_raw_sequence(&raw_data)
        };

        if sequence.len() != length
        {
            return Err(anyhow::anyhow!(
                "Sequence length mismatch for '{}': expected {}, got {}",
                chr_name,
                length,
                sequence.len()
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
        gzi_index: &SimpleGziIndex,
        expected_length: usize,
    ) -> Result<Vec<u8>>
    {
        // Calculate the byte range we need in uncompressed space
        let start_uncompressed = fai_entry.offset;
        let end_uncompressed = fai_entry.offset_for_position(fai_entry.length as u64);

        // Find which GZI blocks overlap with our range
        let entries = gzi_index.entries();

        let start_block_idx = entries
            .partition_point(|&(_, unc)| unc < start_uncompressed)
            .saturating_sub(1);

        let end_block_idx = entries
            .partition_point(|&(_, unc)| unc < end_uncompressed)
            .min(entries.len());

        if start_block_idx >= entries.len()
        {
            return Err(anyhow::anyhow!("No GZI blocks found for sequence"));
        }

        // Fetch all required blocks
        let first_block = entries[start_block_idx];
        let last_block = entries[end_block_idx.saturating_sub(1)];

        let fetch_start = first_block.0;
        let fetch_end = last_block.0 + (64 * 1024); // Estimate end

        let compressed_data = Self::fetch_url_range(url, fetch_start, fetch_end).await?;

        let mut result = Vec::with_capacity(expected_length);

        for block_idx in start_block_idx..end_block_idx
        {
            let (comp_off, unc_off) = entries[block_idx];

            let block_start_in_fetch = (comp_off - fetch_start) as usize;
            let block_end_in_fetch = if block_idx + 1 < entries.len()
            {
                (entries[block_idx + 1].0 - fetch_start) as usize
            }
            else
            {
                compressed_data.len()
            };

            if block_end_in_fetch > compressed_data.len()
            {
                return Err(anyhow::anyhow!("Block extends beyond fetched data"));
            }

            let block_data = &compressed_data[block_start_in_fetch..block_end_in_fetch];
            let decompressed = Self::decompress_bgzf_block(block_data)?;

            let block_start_offset = if unc_off >= start_uncompressed
            {
                0
            }
            else
            {
                (start_uncompressed - unc_off) as usize
            };

            let block_end_offset = if unc_off + decompressed.len() as u64 <= end_uncompressed
            {
                decompressed.len()
            }
            else
            {
                (end_uncompressed - unc_off) as usize
            };

            if block_start_offset < decompressed.len() && block_end_offset > block_start_offset
            {
                let extract_end = block_end_offset.min(decompressed.len());
                result.extend_from_slice(&decompressed[block_start_offset..extract_end]);
            }
        }

        Ok(Self::process_raw_sequence(&result))
    }

    /// Decompress a single BGZF block
    #[cfg(target_arch = "wasm32")]
    fn decompress_bgzf_block(data: &[u8]) -> Result<Vec<u8>>
    {
        use flate2::read::GzDecoder;

        if data.len() < 18
        {
            return Err(anyhow::anyhow!("BGZF block too small: {}", data.len()));
        }

        if data[0] != 0x1f || data[1] != 0x8b
        {
            return Err(anyhow::anyhow!("Invalid gzip header"));
        }

        let mut decoder = GzDecoder::new(data);
        let mut decompressed = Vec::new();
        std::io::Read::read_to_end(&mut decoder, &mut decompressed)
            .map_err(|e| anyhow::anyhow!("Decompression failed: {}", e))?;

        Ok(decompressed)
    }

    /// Process raw sequence bytes (removes newlines)
    #[cfg(target_arch = "wasm32")]
    fn process_raw_sequence(data: &[u8]) -> Vec<u8>
    {
        let mut result = Vec::with_capacity(data.len());
        for &byte in data
        {
            if byte != b'\r' && byte != b'\n' && !byte.is_ascii_whitespace()
            {
                result.push(byte.to_ascii_uppercase());
            }
        }
        result
    }

    /// Load genome from a local file path
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(path: &str) -> Result<Self>
    {
        if path.starts_with("http://") || path.starts_with("https://")
        {
            return Self::from_url(path);
        }

        let path_obj = Path::new(path);

        if let Ok(genome) = Self::try_indexed_local(path_obj)
        {
            return Ok(genome);
        }

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
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_url(url: &str) -> Result<Self>
    {
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
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_url_unindexed(url: &str, block_size: Option<u64>) -> Result<Self>
    {
        let remote_reader = match block_size
        {
            Some(size) => RemoteReader::new(url)
                .map(|r| r.with_block_size(size))
                .with_context(|| format!("Failed to create remote reader for: {}", url))?,
            None => RemoteReader::new(url)
                .with_context(|| format!("Failed to create remote reader for: {}", url))?,
        };

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
            let sequence: Vec<u8> = record
                .seq()
                .iter()
                .map(|&b| b.to_ascii_uppercase())
                .collect();
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

    pub fn get_chromosome_info(&self, name: &str) -> Option<&ChromosomeInfo>
    {
        self.chromosomes.get(name)
    }

    pub fn is_indexed(&self) -> bool
    {
        !matches!(self.data, GenomeData::Full(_))
    }

    pub fn get_chromosome(&self, chr_name: &str) -> Option<&Chromosome>
    {
        match &self.data
        {
            GenomeData::Full(chromosomes) => chromosomes.get(chr_name),
            _ => self.chromosome_cache.get(chr_name),
        }
    }

    pub fn ensure_chromosome_loaded(&mut self, chr_name: &str) -> Result<()>
    {
        match &mut self.data
        {
            GenomeData::Full(_) => Ok(()),
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
                    let sequence: Vec<u8> = record
                        .seq()
                        .iter()
                        .map(|&b| b.to_ascii_uppercase())
                        .collect();

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
                    let sequence: Vec<u8> = record
                        .seq()
                        .iter()
                        .map(|&b| b.to_ascii_uppercase())
                        .collect();

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
                if self.chromosome_cache.contains_key(chr_name)
                {
                    Ok(())
                }
                else
                {
                    Err(anyhow::anyhow!(
                        "WASM indexed genome requires async loading. Use load_chromosome_async() instead."
                    ))
                }
            }
        }
    }

    pub fn get_full_sequence(&mut self, chr_name: &str) -> Result<&[u8]>
    {
        self.ensure_chromosome_loaded(chr_name)?;
        match &self.data
        {
            GenomeData::Full(chromosomes) => chromosomes
                .get(chr_name)
                .map(|c| c.sequence.as_slice())
                .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not found", chr_name)),
            _ => self
                .chromosome_cache
                .get(chr_name)
                .map(|c| c.sequence.as_slice())
                .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not cached", chr_name)),
        }
    }

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
                if let Some(cached) = self.chromosome_cache.get(chr_name)
                {
                    let start = start.min(cached.length);
                    let end = end.min(cached.length);
                    return Ok(cached.sequence[start..end].to_vec());
                }

                let data = reader
                    .fetch_range(chr_name, start as u64, end as u64)
                    .with_context(|| {
                        format!("Failed to fetch range {}:{}-{}", chr_name, start, end)
                    })?;
                Ok(data.iter().map(|&b| b.to_ascii_uppercase()).collect())
            }
            #[cfg(not(target_arch = "wasm32"))]
            GenomeData::IndexedRemote(reader) =>
            {
                if let Some(cached) = self.chromosome_cache.get(chr_name)
                {
                    let start = start.min(cached.length);
                    let end = end.min(cached.length);
                    return Ok(cached.sequence[start..end].to_vec());
                }

                let data = reader
                    .fetch_range(chr_name, start as u64, end as u64)
                    .with_context(|| {
                        format!("Failed to fetch range {}:{}-{}", chr_name, start, end)
                    })?;
                Ok(data.iter().map(|&b| b.to_ascii_uppercase()).collect())
            }
            #[cfg(target_arch = "wasm32")]
            GenomeData::WasmRemoteIndexed(_) =>
            {
                if let Some(cached) = self.chromosome_cache.get(chr_name)
                {
                    let start = start.min(cached.length);
                    let end = end.min(cached.length);
                    return Ok(cached.sequence[start..end].to_vec());
                }
                Err(anyhow::anyhow!(
                    "Chromosome '{}' not loaded. Call load_chromosome_async() first.",
                    chr_name
                ))
            }
        }
    }

    pub fn clear_cache(&mut self)
    {
        self.chromosome_cache.clear();
    }

    pub fn get_sequence(&self, chr: &str, start: usize, end: usize) -> Option<&[u8]>
    {
        match &self.data
        {
            GenomeData::Full(chromosomes) => chromosomes
                .get(chr)
                .and_then(|c| c.sequence.get(start..end)),
            _ => self
                .chromosome_cache
                .get(chr)
                .and_then(|c| c.sequence.get(start..end)),
        }
    }

    pub fn has_sequence(&self, seq_id: &str) -> bool
    {
        self.chromosomes.contains_key(seq_id)
    }

    pub fn sequence_names(&self) -> Vec<String>
    {
        self.chromosomes.keys().cloned().collect()
    }

    pub fn sequence_count(&self) -> usize
    {
        self.chromosomes.len()
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

        gc_count as f32 / region_length(start, end) as f32
    }

    pub fn get_base_at(&self, pos: usize) -> Option<char>
    {
        self.sequence.get(pos).map(|&b| b as char)
    }
}

fn region_length(start: usize, end: usize) -> usize
{
    if end > start { end - start } else { 0 }
}
