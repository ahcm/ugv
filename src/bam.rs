use crate::gff::Strand;
use anyhow::Result;
use std::collections::HashMap;

use noodles::bam;
use noodles::sam;

#[cfg(target_arch = "wasm32")]
use std::io::Cursor;

#[cfg(not(target_arch = "wasm32"))]
use anyhow::Context;
#[cfg(not(target_arch = "wasm32"))]
use std::io::{self, Cursor, Read, Seek, SeekFrom};
#[cfg(not(target_arch = "wasm32"))]
use std::path::{Path, PathBuf};
#[cfg(not(target_arch = "wasm32"))]
use std::sync::Arc;

// Memory management constants
const MAX_LOADED_REGION_SIZE: usize = 10_000_000; // 10Mb max region
const MAX_LOADED_ALIGNMENTS: usize = 100_000_000; // Max 100m reads
const MAX_LOADED_ALIGNMENT_BYTES: usize = 4 * 1024 * 1024 * 1024; // 4 GiB soft cap
const VIEWPORT_BUFFER_MIN: usize = 100_000; // 100Kb min buffer around viewport
const VIEWPORT_BUFFER_MAX: usize = 2_000_000; // 2Mb max buffer around viewport

/// CIGAR operation types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CigarOp
{
    Match(usize),       // M: alignment match (can be match or mismatch)
    Insertion(usize),   // I: insertion to reference
    Deletion(usize),    // D: deletion from reference
    Skip(usize),        // N: skipped region (intron)
    SoftClip(usize),    // S: soft clipping
    HardClip(usize),    // H: hard clipping
    Padding(usize),     // P: padding
    SeqMatch(usize),    // =: sequence match
    SeqMismatch(usize), // X: sequence mismatch
}

/// Variant type
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VariantType
{
    SNP,
    Insertion,
    Deletion,
}

/// Variant detected from alignment
#[derive(Debug, Clone)]
pub struct Variant
{
    pub position: usize,
    pub variant_type: VariantType,
    pub ref_base: Option<u8>,
    pub alt_base: Option<u8>,
    pub quality: u8,
}

/// Methylation modification type
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ModificationType
{
    FiveMC,  // 5-methylcytosine (5mC) - most common
    FiveHMC, // 5-hydroxymethylcytosine (5hmC)
    SixMA,   // 6-methyladenine (6mA)
    Other,
}

/// Single methylation call from a read
#[derive(Debug, Clone)]
pub struct MethylationCall
{
    pub position: usize,      // Genomic position (0-based)
    pub read_position: usize, // Position in read sequence
    pub modification: ModificationType,
    pub probability: u8, // 0-255 (ML tag value)
    pub strand: Strand,
}

/// Aggregated methylation data at a genomic position
#[derive(Debug, Clone, Default)]
pub struct MethylationPoint
{
    pub position: usize,
    pub total_calls: u32,      // Number of reads covering this position
    pub methylated_calls: u32, // Reads with methylation probability > 128
    pub probability_sum: u64,  // Sum of probabilities for averaging
}

/// Single alignment record with pre-computed fields
#[derive(Debug, Clone)]
pub struct AlignmentRecord
{
    pub name: String,
    pub reference_start: usize,            // 0-based position
    pub reference_end: usize,              // Computed from alignment length
    pub sequence: Vec<u8>,                 // Read sequence
    pub quality_scores: Vec<u8>,           // Phred quality scores
    pub mapping_quality: u8,               // MAPQ
    pub flags: u16,                        // SAM flags
    pub strand: Strand,                    // Forward or Reverse
    pub cigar: Vec<CigarOp>,               // Parsed CIGAR operations
    pub variants: Vec<Variant>,            // Pre-computed from CIGAR + MD tag
    pub methylation: Vec<MethylationCall>, // Methylation calls from MM/ML tags
}

impl AlignmentRecord
{
    /// Estimate in-memory footprint (struct + owned heap allocations).
    fn estimated_bytes(&self) -> usize
    {
        std::mem::size_of::<Self>()
            + self.name.capacity()
            + self.sequence.capacity() * std::mem::size_of::<u8>()
            + self.quality_scores.capacity() * std::mem::size_of::<u8>()
            + self.cigar.capacity() * std::mem::size_of::<CigarOp>()
            + self.variants.capacity() * std::mem::size_of::<Variant>()
            + self.methylation.capacity() * std::mem::size_of::<MethylationCall>()
    }

    /// Check if read is mapped
    pub fn is_mapped(&self) -> bool
    {
        (self.flags & 0x4) == 0
    }

    /// Check if read is on reverse strand
    pub fn is_reverse(&self) -> bool
    {
        (self.flags & 0x10) != 0
    }

    /// Extract variants from CIGAR string
    pub fn extract_variants(&mut self)
    {
        let mut variants = Vec::new();
        let mut ref_pos = self.reference_start;
        let mut read_pos = 0;

        for cigar_op in &self.cigar
        {
            match cigar_op
            {
                CigarOp::SeqMismatch(len) =>
                {
                    // Explicit mismatches from = and X CIGAR ops
                    for i in 0..*len
                    {
                        if read_pos + i < self.sequence.len()
                        {
                            variants.push(Variant {
                                position: ref_pos + i,
                                variant_type: VariantType::SNP,
                                ref_base: None,
                                alt_base: Some(self.sequence[read_pos + i]),
                                quality: self
                                    .quality_scores
                                    .get(read_pos + i)
                                    .copied()
                                    .unwrap_or(0),
                            });
                        }
                    }
                    ref_pos += len;
                    read_pos += len;
                }
                CigarOp::Insertion(len) =>
                {
                    variants.push(Variant {
                        position: ref_pos,
                        variant_type: VariantType::Insertion,
                        ref_base: None,
                        alt_base: self.sequence.get(read_pos).copied(),
                        quality: self.quality_scores.get(read_pos).copied().unwrap_or(0),
                    });
                    read_pos += len;
                }
                CigarOp::Deletion(len) =>
                {
                    variants.push(Variant {
                        position: ref_pos,
                        variant_type: VariantType::Deletion,
                        ref_base: None,
                        alt_base: None,
                        quality: 0,
                    });
                    ref_pos += len;
                }
                CigarOp::Match(len) | CigarOp::SeqMatch(len) =>
                {
                    ref_pos += len;
                    read_pos += len;
                }
                CigarOp::SoftClip(len) =>
                {
                    read_pos += len;
                }
                CigarOp::Skip(len) | CigarOp::HardClip(len) | CigarOp::Padding(len) =>
                {
                    ref_pos += len;
                }
            }
        }

        self.variants = variants;
    }

    /// Compute reference length from CIGAR operations
    fn compute_reference_length(cigar_ops: &[CigarOp]) -> usize
    {
        let mut length = 0;
        for op in cigar_ops
        {
            match op
            {
                CigarOp::Match(len)
                | CigarOp::Deletion(len)
                | CigarOp::Skip(len)
                | CigarOp::SeqMatch(len)
                | CigarOp::SeqMismatch(len) =>
                {
                    length += len;
                }
                _ =>
                {}
            }
        }
        length
    }
}

/// Coverage at specific position
#[derive(Debug, Clone, Default)]
pub struct CoveragePoint
{
    pub position: usize,
    pub depth: u32,
    pub forward_depth: u32,
    pub reverse_depth: u32,
    pub quality_sum: u64,
}

/// Track for one reference/chromosome
#[derive(Clone)]
pub struct AlignmentTrack
{
    pub reference_name: String,
    pub reference_length: usize,
    pub records: Vec<AlignmentRecord>,
    pub coverage: Vec<CoveragePoint>,
    pub methylation: Vec<MethylationPoint>,
    pub loaded_region: Option<(usize, usize)>,
}

impl AlignmentTrack
{
    pub fn new(reference_name: String, reference_length: usize) -> Self
    {
        Self {
            reference_name,
            reference_length,
            records: Vec::new(),
            coverage: Vec::new(),
            methylation: Vec::new(),
            loaded_region: None,
        }
    }

    /// Compute coverage histogram with binning
    pub fn compute_coverage(&mut self, bin_size: usize)
    {
        if self.reference_length == 0 || self.loaded_region.is_none()
        {
            return;
        }

        let (region_start, region_end) = self.loaded_region.unwrap();
        if region_end <= region_start
        {
            return;
        }
        let region_length = region_end - region_start;
        let num_bins = (region_length + bin_size - 1) / bin_size;
        let mut coverage = vec![CoveragePoint::default(); num_bins];

        for (bin_idx, point) in coverage.iter_mut().enumerate()
        {
            point.position = region_start + bin_idx * bin_size;
        }

        for record in &self.records
        {
            if record.reference_start < region_end && record.reference_end > region_start
            {
                let start_pos = record.reference_start.max(region_start) - region_start;
                let end_pos = record.reference_end.min(region_end) - region_start;

                let start_bin = start_pos / bin_size;
                let end_bin = end_pos / bin_size;

                for bin_idx in start_bin..=end_bin.min(num_bins - 1)
                {
                    coverage[bin_idx].depth += 1;
                    coverage[bin_idx].quality_sum += record.mapping_quality as u64;

                    match record.strand
                    {
                        Strand::Forward => coverage[bin_idx].forward_depth += 1,
                        Strand::Reverse => coverage[bin_idx].reverse_depth += 1,
                        _ =>
                        {}
                    }
                }
            }
        }

        self.coverage = coverage;
    }

    /// Compute methylation levels with binning
    pub fn compute_methylation(&mut self, bin_size: usize)
    {
        if self.reference_length == 0 || self.loaded_region.is_none()
        {
            return;
        }

        let (region_start, region_end) = self.loaded_region.unwrap();
        let region_length = region_end - region_start;
        let num_bins = (region_length + bin_size - 1) / bin_size;
        let mut methylation = vec![MethylationPoint::default(); num_bins];

        // Initialize positions
        for (bin_idx, point) in methylation.iter_mut().enumerate()
        {
            point.position = region_start + bin_idx * bin_size;
        }

        // Aggregate methylation calls from all reads
        for record in &self.records
        {
            for call in &record.methylation
            {
                if call.position >= region_start && call.position < region_end
                {
                    let bin_idx = (call.position - region_start) / bin_size;
                    if bin_idx < num_bins
                    {
                        methylation[bin_idx].total_calls += 1;
                        methylation[bin_idx].probability_sum += call.probability as u64;
                        if call.probability > 128
                        {
                            methylation[bin_idx].methylated_calls += 1;
                        }
                    }
                }
            }
        }

        self.methylation = methylation;
    }

    /// Clear loaded data to free memory
    pub fn clear(&mut self)
    {
        self.records.clear();
        self.coverage.clear();
        self.methylation.clear();
        self.loaded_region = None;
    }
}

/// Container for BAM data with lazy loading support
#[derive(Clone)]
pub struct AlignmentData
{
    #[cfg(not(target_arch = "wasm32"))]
    pub bam_path: String,
    #[cfg(not(target_arch = "wasm32"))]
    _remote_storage: Option<Arc<RemoteBamStorage>>,
    #[cfg(not(target_arch = "wasm32"))]
    remote_source: Option<Arc<RemoteBamSource>>,

    #[cfg(target_arch = "wasm32")]
    pub bam_bytes: Vec<u8>,

    pub header: sam::Header,
    pub reference_lengths: HashMap<String, usize>,

    // Currently loaded tracks (only for displayed region)
    pub loaded_tracks: HashMap<String, AlignmentTrack>,
    pub previous_tracks: HashMap<String, AlignmentTrack>,
}

#[cfg(not(target_arch = "wasm32"))]
#[derive(Debug)]
struct RemoteBamStorage
{
    dir: PathBuf,
}

#[cfg(not(target_arch = "wasm32"))]
#[derive(Debug, Clone)]
struct RemoteBamSource
{
    url: String,
    bai_data: Vec<u8>,
    content_length: u64,
    client: reqwest::blocking::Client,
}

#[cfg(not(target_arch = "wasm32"))]
struct HttpRangeReader
{
    client: reqwest::blocking::Client,
    url: String,
    position: u64,
    content_length: u64,
    cache_start: u64,
    cache_data: Vec<u8>,
}

#[cfg(not(target_arch = "wasm32"))]
impl HttpRangeReader
{
    const FETCH_WINDOW_BYTES: u64 = 30_048_576; // 1 MiB

    fn new(client: reqwest::blocking::Client, url: String, content_length: u64) -> Self
    {
        Self {
            client,
            url,
            position: 0,
            content_length,
            cache_start: 0,
            cache_data: Vec::new(),
        }
    }

    fn cache_end(&self) -> u64
    {
        self.cache_start + self.cache_data.len() as u64
    }

    fn refill_cache(&mut self, min_bytes: usize) -> io::Result<()>
    {
        if self.position >= self.content_length
        {
            self.cache_data.clear();
            return Ok(());
        }

        let want = std::cmp::max(min_bytes as u64, Self::FETCH_WINDOW_BYTES);
        let start = self.position;
        let end = (start + want - 1).min(self.content_length - 1);
        let range_value = format!("bytes={}-{}", start, end);

        let response = self
            .client
            .get(&self.url)
            .header(reqwest::header::RANGE, range_value)
            .send()
            .map_err(io::Error::other)?;

        let status = response.status();
        if !(status == reqwest::StatusCode::PARTIAL_CONTENT || status == reqwest::StatusCode::OK)
        {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("unexpected HTTP status {} while reading {}", status, self.url),
            ));
        }

        let bytes = response.bytes().map_err(io::Error::other)?;
        self.cache_start = start;
        self.cache_data = bytes.to_vec();
        Ok(())
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl Read for HttpRangeReader
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize>
    {
        if buf.is_empty()
        {
            return Ok(0);
        }
        if self.position >= self.content_length
        {
            return Ok(0);
        }

        if self.position < self.cache_start || self.position >= self.cache_end()
        {
            self.refill_cache(buf.len())?;
        }

        if self.cache_data.is_empty()
        {
            return Ok(0);
        }

        let cache_offset = (self.position - self.cache_start) as usize;
        if cache_offset >= self.cache_data.len()
        {
            self.refill_cache(buf.len())?;
            if self.cache_data.is_empty()
            {
                return Ok(0);
            }
        }

        let available = &self.cache_data[(self.position - self.cache_start) as usize..];
        let to_copy = std::cmp::min(buf.len(), available.len());
        buf[..to_copy].copy_from_slice(&available[..to_copy]);
        self.position += to_copy as u64;
        Ok(to_copy)
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl Seek for HttpRangeReader
{
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64>
    {
        let new_pos: i128 = match pos
        {
            SeekFrom::Start(offset) => offset as i128,
            SeekFrom::Current(delta) => self.position as i128 + delta as i128,
            SeekFrom::End(delta) => self.content_length as i128 + delta as i128,
        };

        if new_pos < 0
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid seek to negative position",
            ));
        }

        self.position = (new_pos as u64).min(self.content_length);
        Ok(self.position)
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl Drop for RemoteBamStorage
{
    fn drop(&mut self)
    {
        if let Err(e) = std::fs::remove_dir_all(&self.dir)
        {
            eprintln!(
                "Warning: failed to remove temporary BAM directory {}: {}",
                self.dir.display(),
                e
            );
        }
    }
}

impl AlignmentData
{
    fn hit_alignment_limit(count: usize, loaded_bytes: usize, next_alignment_bytes: usize) -> bool
    {
        count >= MAX_LOADED_ALIGNMENTS
            || loaded_bytes.saturating_add(next_alignment_bytes) > MAX_LOADED_ALIGNMENT_BYTES
    }

    fn extract_reference_lengths(header: &sam::Header) -> HashMap<String, usize>
    {
        header
            .reference_sequences()
            .iter()
            .map(|(name, reference_sequence)| {
                (name.to_string(), usize::from(reference_sequence.length()))
            })
            .collect()
    }

    fn reference_length(&self, chromosome: &str) -> usize
    {
        *self.reference_lengths.get(chromosome).unwrap_or(&0)
    }

    /// Create from file path (native only) - does NOT load records yet
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(bam_path: &str) -> Result<Self>
    {
        use std::fs::File;

        let mut reader = bam::io::Reader::new(File::open(bam_path)?);
        let header = reader.read_header()?;
        let reference_lengths = Self::extract_reference_lengths(&header);

        Ok(AlignmentData {
            bam_path: bam_path.to_string(),
            _remote_storage: None,
            remote_source: None,
            header,
            reference_lengths,
            loaded_tracks: HashMap::new(),
            previous_tracks: HashMap::new(),
        })
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_path_or_url(path_or_url: &str) -> Result<Self>
    {
        if path_or_url.starts_with("http://") || path_or_url.starts_with("https://")
        {
            return Self::from_url(path_or_url);
        }
        Self::from_file(path_or_url)
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn from_url(url: &str) -> Result<Self>
    {
        use reqwest::{StatusCode, header};
        use std::fs::File;
        use std::io::copy;
        use std::time::{SystemTime, UNIX_EPOCH};

        fn fetch_content_length(client: &reqwest::blocking::Client, source_url: &str)
        -> Result<u64>
        {
            // Prefer HEAD, but fall back to GET range probe when servers do not support HEAD.
            if let Ok(response) = client.head(source_url).send()
            {
                if response.status().is_success()
                {
                    if let Some(len) = response
                        .headers()
                        .get(header::CONTENT_LENGTH)
                        .and_then(|v| v.to_str().ok())
                        .and_then(|v| v.parse::<u64>().ok())
                    {
                        return Ok(len);
                    }
                }
            }

            let response = client
                .get(source_url)
                .header(header::RANGE, "bytes=0-0")
                .send()
                .with_context(|| format!("Failed to probe content length for {}", source_url))?;

            if !response.status().is_success() && response.status() != StatusCode::PARTIAL_CONTENT
            {
                return Err(anyhow::anyhow!(
                    "HTTP error {} while probing content length for {}",
                    response.status(),
                    source_url
                ));
            }

            if let Some(content_range) = response
                .headers()
                .get(header::CONTENT_RANGE)
                .and_then(|v| v.to_str().ok())
            {
                // Format: bytes 0-0/12345
                if let Some(total) = content_range.split('/').nth(1)
                {
                    if let Ok(len) = total.parse::<u64>()
                    {
                        return Ok(len);
                    }
                }
            }

            if let Some(len) = response
                .headers()
                .get(header::CONTENT_LENGTH)
                .and_then(|v| v.to_str().ok())
                .and_then(|v| v.parse::<u64>().ok())
            {
                return Ok(len);
            }

            Err(anyhow::anyhow!("Could not determine content length for {}", source_url))
        }

        fn download_to_path(
            client: &reqwest::blocking::Client,
            source_url: &str,
            dest: &Path,
        ) -> Result<()>
        {
            let mut response = client
                .get(source_url)
                .send()
                .with_context(|| format!("Failed to download URL: {}", source_url))?;

            if !response.status().is_success()
            {
                return Err(anyhow::anyhow!(
                    "HTTP error {} while downloading {}",
                    response.status(),
                    source_url
                ));
            }

            let mut file = File::create(dest)
                .with_context(|| format!("Failed to create file {}", dest.display()))?;
            copy(&mut response, &mut file)
                .with_context(|| format!("Failed to write downloaded file {}", dest.display()))?;
            Ok(())
        }

        fn download_optional_bai_bytes(
            client: &reqwest::blocking::Client,
            source_url: &str,
        ) -> Result<Option<Vec<u8>>>
        {
            let response = client
                .get(source_url)
                .send()
                .with_context(|| format!("Failed to check BAM index URL: {}", source_url))?;

            if response.status() == StatusCode::NOT_FOUND
            {
                return Ok(None);
            }
            if !response.status().is_success()
            {
                return Ok(None);
            }

            let bytes = response.bytes().with_context(|| {
                format!("Failed to read BAM index response body: {}", source_url)
            })?;
            Ok(Some(bytes.to_vec()))
        }

        let client = reqwest::blocking::Client::builder()
            .build()
            .context("Failed to create HTTP client for BAM download")?;

        let mut bai_candidates = vec![format!("{}.bai", url)];
        if let Some(stripped) = url.strip_suffix(".bam")
        {
            bai_candidates.push(format!("{}.bai", stripped));
        }

        let mut bai_data = None;
        for candidate in bai_candidates
        {
            if let Some(index_bytes) = download_optional_bai_bytes(&client, &candidate)?
            {
                bai_data = Some(index_bytes);
                break;
            }
        }

        let has_remote_index = bai_data.is_some();
        if let Some(bai_data) = bai_data
        {
            let ranged_result: Result<Self> = (|| {
                let content_length = fetch_content_length(&client, url)?;

                let mut index_reader = bam::bai::io::Reader::new(Cursor::new(&bai_data));
                let index = index_reader.read_index()?;

                let range_reader =
                    HttpRangeReader::new(client.clone(), url.to_string(), content_length);
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .set_index(index)
                    .build_from_reader(range_reader)?;

                let header = reader.read_header()?;
                let reference_lengths = Self::extract_reference_lengths(&header);

                Ok(AlignmentData {
                    bam_path: url.to_string(),
                    _remote_storage: None,
                    remote_source: Some(Arc::new(RemoteBamSource {
                        url: url.to_string(),
                        bai_data: bai_data.clone(),
                        content_length,
                        client: client.clone(),
                    })),
                    header,
                    reference_lengths,
                    loaded_tracks: HashMap::new(),
                    previous_tracks: HashMap::new(),
                })
            })();

            match ranged_result
            {
                Ok(data) => return Ok(data),
                Err(e) =>
                {
                    eprintln!(
                        "Warning: failed to initialize ranged BAM reader for {} ({}); falling back to full download",
                        url, e
                    );
                }
            }
        }

        if !has_remote_index
        {
            eprintln!(
                "Warning: no remote BAM index found for {} (.bam.bai/.bai); falling back to full download",
                url
            );
        }

        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let temp_dir =
            std::env::temp_dir().join(format!("ugv-bam-{}-{}", std::process::id(), unique));
        std::fs::create_dir_all(&temp_dir)
            .with_context(|| format!("Failed to create temp directory {}", temp_dir.display()))?;

        let storage = Arc::new(RemoteBamStorage {
            dir: temp_dir.clone(),
        });
        let bam_path = temp_dir.join("remote.bam");
        download_to_path(&client, url, &bam_path)?;

        // Keep compatibility with build_from_path index auto-discovery if a sidecar index appears.
        if let Some(index_bytes) = download_optional_bai_bytes(&client, &format!("{}.bai", url))?
        {
            let bam_dot_bai = PathBuf::from(format!("{}.bai", bam_path.display()));
            let bai_ext = bam_path.with_extension("bai");
            let _ = std::fs::write(&bam_dot_bai, &index_bytes);
            let _ = std::fs::write(&bai_ext, &index_bytes);
        }

        let mut reader = bam::io::Reader::new(File::open(&bam_path)?);
        let header = reader.read_header()?;
        let reference_lengths = Self::extract_reference_lengths(&header);

        Ok(AlignmentData {
            bam_path: bam_path.display().to_string(),
            _remote_storage: Some(storage),
            remote_source: None,
            header,
            reference_lengths,
            loaded_tracks: HashMap::new(),
            previous_tracks: HashMap::new(),
        })
    }

    /// Create from bytes (WASM) - does NOT load records yet
    #[cfg(target_arch = "wasm32")]
    pub fn from_bytes(data: Vec<u8>) -> Result<Self>
    {
        let cursor = Cursor::new(&data);
        let mut reader = bam::io::Reader::new(cursor);
        let header = reader.read_header()?;
        let reference_lengths = Self::extract_reference_lengths(&header);

        Ok(AlignmentData {
            bam_bytes: data,
            header,
            reference_lengths,
            loaded_tracks: HashMap::new(),
            previous_tracks: HashMap::new(),
        })
    }

    fn buffered_region(
        &self,
        viewport_start: usize,
        viewport_end: usize,
        reference_length: usize,
    ) -> (usize, usize)
    {
        let viewport_size = viewport_end.saturating_sub(viewport_start);
        let buffer = (viewport_size / 2).clamp(VIEWPORT_BUFFER_MIN, VIEWPORT_BUFFER_MAX);
        let buffer_start = viewport_start.saturating_sub(buffer);
        let buffer_end = if reference_length == 0
        {
            viewport_end.saturating_add(buffer)
        }
        else
        {
            (viewport_end + buffer).min(reference_length)
        };
        (buffer_start, buffer_end)
    }

    /// Query a specific region - loads from BAM if not already cached
    #[cfg(not(target_arch = "wasm32"))]
    pub fn query_region(
        &mut self,
        chromosome: &str,
        viewport_start: usize,
        viewport_end: usize,
    ) -> Result<Option<&AlignmentTrack>>
    {
        let reference_length = self.reference_length(chromosome);
        let (buffer_start, buffer_end) =
            self.buffered_region(viewport_start, viewport_end, reference_length);

        // Limit region size
        let region_size = buffer_end - buffer_start;
        if region_size > MAX_LOADED_REGION_SIZE
        {
            // Region too large, load only viewport without buffer
            return self.query_region_exact(chromosome, viewport_start, viewport_end);
        }

        // Check if we need to load (avoid borrow checker issues)
        let mut need_load = true;
        if let Some(track) = self.loaded_tracks.get(chromosome)
        {
            if let Some((loaded_start, loaded_end)) = track.loaded_region
            {
                // Check if requested region is outside loaded region
                need_load = buffer_start < loaded_start || buffer_end > loaded_end;
            }
        }
        else if let Some(track) = self.previous_tracks.get(chromosome)
        {
            if let Some((loaded_start, loaded_end)) = track.loaded_region
            {
                if buffer_start >= loaded_start && buffer_end <= loaded_end
                {
                    if let Some(track) = self.previous_tracks.remove(chromosome)
                    {
                        self.loaded_tracks.insert(chromosome.to_string(), track);
                        need_load = false;
                    }
                }
            }
        };

        // Load if needed
        if need_load
        {
            self.load_region(chromosome, buffer_start, buffer_end)?;
        }

        // Now return the reference
        Ok(self.loaded_tracks.get(chromosome))
    }

    /// Query exact region without buffer (for large viewports)
    #[cfg(not(target_arch = "wasm32"))]
    fn query_region_exact(
        &mut self,
        chromosome: &str,
        start: usize,
        end: usize,
    ) -> Result<Option<&AlignmentTrack>>
    {
        self.load_region(chromosome, start, end)?;
        Ok(self.loaded_tracks.get(chromosome))
    }

    /// Load a specific region from BAM file using indexed reader
    #[cfg(not(target_arch = "wasm32"))]
    fn load_region(&mut self, chromosome: &str, start: usize, end: usize) -> Result<()>
    {
        use noodles::core::Region;

        let reference_length = *self.reference_lengths.get(chromosome).unwrap_or(&0);
        let mut track = AlignmentTrack::new(chromosome.to_string(), reference_length);
        track.loaded_region = Some((start, end));

        if let Some(remote_source) = self.remote_source.clone()
        {
            let mut index_reader =
                bam::bai::io::Reader::new(Cursor::new(remote_source.bai_data.as_slice()));
            let index = index_reader.read_index()?;

            let range_reader = HttpRangeReader::new(
                remote_source.client.clone(),
                remote_source.url.clone(),
                remote_source.content_length,
            );
            let mut reader = bam::io::indexed_reader::Builder::default()
                .set_index(index)
                .build_from_reader(range_reader)?;

            let region = format!("{}:{}-{}", chromosome, start + 1, end);
            if let Ok(region) = region.parse::<Region>()
            {
                let header = reader.read_header()?;
                let query = reader.query(&header, &region)?;

                let mut count = 0;
                let mut loaded_bytes = 0usize;
                for result in query.records()
                {
                    let record = result?;
                    if let Some(alignment) = Self::parse_record(&record)?
                    {
                        let alignment_bytes = alignment.estimated_bytes();
                        if Self::hit_alignment_limit(count, loaded_bytes, alignment_bytes)
                        {
                            eprintln!(
                                "Warning: Hit loading limit for region (max {} alignments or {} MiB estimated)",
                                MAX_LOADED_ALIGNMENTS,
                                MAX_LOADED_ALIGNMENT_BYTES / (1024 * 1024)
                            );
                            break;
                        }
                        track.records.push(alignment);
                        count += 1;
                        loaded_bytes += alignment_bytes;
                    }
                }
            }
        }
        else
        {
            // Try local indexed reader first.
            let indexed_result =
                bam::io::indexed_reader::Builder::default().build_from_path(&self.bam_path);

            if let Ok(mut reader) = indexed_result
            {
                let region = format!("{}:{}-{}", chromosome, start + 1, end);
                if let Ok(region) = region.parse::<Region>()
                {
                    let header = reader.read_header()?;
                    let query = reader.query(&header, &region)?;

                    let mut count = 0;
                    let mut loaded_bytes = 0usize;
                    for result in query.records()
                    {
                        let record = result?;
                        if let Some(alignment) = Self::parse_record(&record)?
                        {
                            let alignment_bytes = alignment.estimated_bytes();
                            if Self::hit_alignment_limit(count, loaded_bytes, alignment_bytes)
                            {
                                eprintln!(
                                    "Warning: Hit loading limit for region (max {} alignments or {} MiB estimated)",
                                    MAX_LOADED_ALIGNMENTS,
                                    MAX_LOADED_ALIGNMENT_BYTES / (1024 * 1024)
                                );
                                break;
                            }
                            track.records.push(alignment);
                            count += 1;
                            loaded_bytes += alignment_bytes;
                        }
                    }
                }
            }
            else
            {
                // Fallback: no index available, scan entire file (slow)
                eprintln!(
                    "Warning: BAM index not found for {}, scanning entire file (slow)",
                    self.bam_path
                );

                use std::fs::File;
                let mut reader = bam::io::Reader::new(File::open(&self.bam_path)?);
                let header = reader.read_header()?;

                let mut count = 0;
                let mut loaded_bytes = 0usize;
                for result in reader.records()
                {
                    let record = result?;

                    // Check if record is in our region
                    if let Some(Ok(ref_seq_id)) = record.reference_sequence_id()
                    {
                        if let Some((ref_name, _)) =
                            header.reference_sequences().get_index(ref_seq_id)
                        {
                            if ref_name.to_string() == chromosome
                            {
                                if let Some(Ok(pos)) = record.alignment_start()
                                {
                                    let pos_0based = usize::from(pos) - 1;
                                    if pos_0based >= start && pos_0based <= end
                                    {
                                        if let Some(alignment) = Self::parse_record(&record)?
                                        {
                                            let alignment_bytes = alignment.estimated_bytes();
                                            if Self::hit_alignment_limit(
                                                count,
                                                loaded_bytes,
                                                alignment_bytes,
                                            )
                                            {
                                                eprintln!(
                                                    "Warning: Hit loading limit for region (max {} alignments or {} MiB estimated)",
                                                    MAX_LOADED_ALIGNMENTS,
                                                    MAX_LOADED_ALIGNMENT_BYTES / (1024 * 1024)
                                                );
                                                break;
                                            }
                                            track.records.push(alignment);
                                            count += 1;
                                            loaded_bytes += alignment_bytes;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Compute coverage and methylation for this region
        track.compute_coverage(100); // 100bp bins
        track.compute_methylation(100); // 100bp bins

        // Store in loaded_tracks, replacing any previous data for this chromosome
        if let Some(existing) = self.loaded_tracks.remove(chromosome)
        {
            self.previous_tracks
                .insert(chromosome.to_string(), existing);
        }
        self.loaded_tracks.insert(chromosome.to_string(), track);

        Ok(())
    }

    /// Query region for WASM (no index support, must scan)
    #[cfg(target_arch = "wasm32")]
    pub fn query_region(
        &mut self,
        chromosome: &str,
        _viewport_start: usize,
        _viewport_end: usize,
    ) -> Result<Option<&AlignmentTrack>>
    {
        // For WASM, we'll load once and cache
        // Check if we need to load
        let need_load = if let Some(track) = self.loaded_tracks.get(chromosome)
        {
            track.loaded_region.is_none()
        }
        else
        {
            true
        };

        // Load if needed
        if need_load
        {
            self.load_chromosome_wasm(chromosome)?;
        }

        // Return the reference
        Ok(self.loaded_tracks.get(chromosome))
    }

    /// Load entire chromosome for WASM (no index available in browser)
    #[cfg(target_arch = "wasm32")]
    fn load_chromosome_wasm(&mut self, chromosome: &str) -> Result<()>
    {
        let cursor = Cursor::new(&self.bam_bytes);
        let mut reader = bam::io::Reader::new(cursor);
        let header = reader.read_header()?;

        let reference_length = self.reference_length(chromosome);
        let mut track = AlignmentTrack::new(chromosome.to_string(), reference_length);
        track.loaded_region = Some((0, reference_length));

        let mut count = 0;
        let mut loaded_bytes = 0usize;
        for result in reader.records()
        {
            let record = result?;

            // Check if this record belongs to our chromosome
            if let Some(Ok(ref_seq_id)) = record.reference_sequence_id()
            {
                if let Some((ref_name, _)) = header.reference_sequences().get_index(ref_seq_id)
                {
                    if ref_name.to_string() == chromosome
                    {
                        if let Some(alignment) = Self::parse_record(&record)?
                        {
                            let alignment_bytes = alignment.estimated_bytes();
                            if Self::hit_alignment_limit(count, loaded_bytes, alignment_bytes)
                            {
                                eprintln!(
                                    "Warning: Hit loading limit for chromosome (max {} alignments or {} MiB estimated)",
                                    MAX_LOADED_ALIGNMENTS,
                                    MAX_LOADED_ALIGNMENT_BYTES / (1024 * 1024)
                                );
                                break;
                            }
                            track.records.push(alignment);
                            count += 1;
                            loaded_bytes += alignment_bytes;
                        }
                    }
                }
            }
        }

        track.compute_coverage(100);
        track.compute_methylation(100);
        self.loaded_tracks.insert(chromosome.to_string(), track);

        Ok(())
    }

    /// Parse methylation information from MM and ML tags
    /// MM tag format: "C+m,0,5,3;" means methylated C at deltas 0, 5, 3 from each C
    /// ML tag: array of probabilities (0-255) for each modification
    fn parse_methylation_tags(
        record: &bam::Record,
        reference_start: usize,
        cigar: &[CigarOp],
        sequence: &[u8],
        strand: Strand,
    ) -> Vec<MethylationCall>
    {
        use noodles::sam::alignment::record::data::field::Tag;

        let mut calls = Vec::new();

        // Try to get MM tag (base modifications)
        let mm_data = record.data();
        let mm_value = mm_data.get(&Tag::BASE_MODIFICATIONS);

        // Try to get ML tag (modification probabilities)
        let ml_value = mm_data.get(&Tag::BASE_MODIFICATION_PROBABILITIES);

        // Parse MM tag if present
        let mm_string = match mm_value
        {
            Some(Ok(value)) =>
            {
                // Convert to string representation
                format!("{:?}", value)
            }
            _ => return calls,
        };

        // Parse ML tag probabilities if present
        let probabilities: Vec<u8> = match ml_value
        {
            Some(Ok(value)) =>
            {
                // Try to extract as array of u8
                let ml_str = format!("{:?}", value);
                // Parse array format like "Array(UInt8([128, 200, 150]))"
                Self::parse_ml_array(&ml_str)
            }
            _ => Vec::new(),
        };

        // Parse MM tag format: "C+m,0,5,3;A+a,1,2;"
        // Each section defines: BaseCode+ModCode,delta1,delta2,...;
        let mm_clean = mm_string
            .trim_matches('"')
            .replace("String(\"", "")
            .replace("\")", "");

        let mut prob_idx = 0;

        for section in mm_clean.split(';')
        {
            if section.is_empty()
            {
                continue;
            }

            let parts: Vec<&str> = section.split(',').collect();
            if parts.is_empty()
            {
                continue;
            }

            // Parse base and modification type (e.g., "C+m" or "C+h")
            let (base, mod_type) = Self::parse_mod_code(parts[0]);
            if base == 0
            {
                continue;
            }

            // Find all positions of the target base in the sequence
            let base_positions: Vec<usize> = sequence
                .iter()
                .enumerate()
                .filter(|(_, &b)| b.to_ascii_uppercase() == base)
                .map(|(i, _)| i)
                .collect();

            // Parse delta values and convert to actual positions
            let mut base_idx = 0;
            for delta_str in parts.iter().skip(1)
            {
                if let Ok(delta) = delta_str.parse::<usize>()
                {
                    base_idx += delta;
                    if base_idx < base_positions.len()
                    {
                        let read_pos = base_positions[base_idx];

                        // Convert read position to reference position using CIGAR
                        if let Some(ref_pos) =
                            Self::read_to_ref_position(read_pos, reference_start, cigar)
                        {
                            let probability = if prob_idx < probabilities.len()
                            {
                                probabilities[prob_idx]
                            }
                            else
                            {
                                200 // Default high probability if ML tag missing
                            };

                            calls.push(MethylationCall {
                                position: ref_pos,
                                read_position: read_pos,
                                modification: mod_type,
                                probability,
                                strand,
                            });
                        }
                        prob_idx += 1;
                    }
                    base_idx += 1; // Move past current base
                }
            }
        }

        calls
    }

    /// Parse modification code from MM tag (e.g., "C+m" -> (b'C', ModificationType::FiveMC))
    fn parse_mod_code(code: &str) -> (u8, ModificationType)
    {
        let code = code.trim();
        if code.len() < 3
        {
            return (0, ModificationType::Other);
        }

        let base = code.chars().next().unwrap_or(' ').to_ascii_uppercase() as u8;
        let mod_char = code.chars().last().unwrap_or(' ');

        let mod_type = match (base, mod_char)
        {
            (b'C', 'm') => ModificationType::FiveMC,  // 5mC
            (b'C', 'h') => ModificationType::FiveHMC, // 5hmC
            (b'A', 'a') => ModificationType::SixMA,   // 6mA
            _ => ModificationType::Other,
        };

        (base, mod_type)
    }

    /// Parse ML array from debug string format
    fn parse_ml_array(ml_str: &str) -> Vec<u8>
    {
        // Handle format like "Array(UInt8([128, 200, 150]))" or similar
        let mut result = Vec::new();

        // Find the innermost array content
        if let Some(start) = ml_str.find('[')
        {
            if let Some(end) = ml_str.rfind(']')
            {
                let array_content = &ml_str[start + 1..end];
                for num_str in array_content.split(',')
                {
                    if let Ok(num) = num_str.trim().parse::<u8>()
                    {
                        result.push(num);
                    }
                }
            }
        }

        result
    }

    /// Convert read position to reference position using CIGAR
    fn read_to_ref_position(
        read_pos: usize,
        reference_start: usize,
        cigar: &[CigarOp],
    ) -> Option<usize>
    {
        let mut ref_pos = reference_start;
        let mut read_offset = 0;

        for op in cigar
        {
            match op
            {
                CigarOp::Match(len) | CigarOp::SeqMatch(len) | CigarOp::SeqMismatch(len) =>
                {
                    if read_offset + len > read_pos
                    {
                        return Some(ref_pos + (read_pos - read_offset));
                    }
                    read_offset += len;
                    ref_pos += len;
                }
                CigarOp::Insertion(len) | CigarOp::SoftClip(len) =>
                {
                    if read_offset + len > read_pos
                    {
                        // Position is within insertion - no reference position
                        return None;
                    }
                    read_offset += len;
                }
                CigarOp::Deletion(len) | CigarOp::Skip(len) =>
                {
                    ref_pos += len;
                }
                CigarOp::HardClip(_) | CigarOp::Padding(_) =>
                {
                    // No effect on positions
                }
            }
        }

        None
    }

    /// Parse a single BAM record into AlignmentRecord
    fn parse_record(record: &bam::Record) -> Result<Option<AlignmentRecord>>
    {
        // Skip unmapped reads
        let flags = record.flags();
        if flags.is_unmapped()
        {
            return Ok(None);
        }

        // Get alignment start (convert from 1-based to 0-based)
        let reference_start = match record.alignment_start()
        {
            Some(Ok(pos)) => usize::from(pos) - 1,
            _ => return Ok(None),
        };

        // Parse CIGAR operations
        let mut cigar_ops = Vec::new();
        for cigar_result in record.cigar().iter()
        {
            let op = cigar_result?;
            use sam::alignment::record::cigar::op::Kind;

            let cigar_op = match op.kind()
            {
                Kind::Match => CigarOp::Match(op.len()),
                Kind::Insertion => CigarOp::Insertion(op.len()),
                Kind::Deletion => CigarOp::Deletion(op.len()),
                Kind::Skip => CigarOp::Skip(op.len()),
                Kind::SoftClip => CigarOp::SoftClip(op.len()),
                Kind::HardClip => CigarOp::HardClip(op.len()),
                Kind::Pad => CigarOp::Padding(op.len()),
                Kind::SequenceMatch => CigarOp::SeqMatch(op.len()),
                Kind::SequenceMismatch => CigarOp::SeqMismatch(op.len()),
            };
            cigar_ops.push(cigar_op);
        }

        // Compute reference end from CIGAR
        let ref_length = AlignmentRecord::compute_reference_length(&cigar_ops);
        let reference_end = reference_start + ref_length;

        // Extract sequence
        let sequence = record.sequence().as_ref().to_vec();

        // Extract quality scores
        let quality_scores: Vec<u8> = record.quality_scores().as_ref().to_vec();

        // Get mapping quality
        let mapping_quality = record.mapping_quality().map(|mq| u8::from(mq)).unwrap_or(0);

        // Get read name
        let name = record
            .name()
            .map(|n| n.to_string())
            .unwrap_or_else(|| "unknown".to_string());

        // Determine strand from flags
        let strand = if flags.is_reverse_complemented()
        {
            Strand::Reverse
        }
        else
        {
            Strand::Forward
        };

        // Parse methylation from MM/ML tags
        let methylation =
            Self::parse_methylation_tags(record, reference_start, &cigar_ops, &sequence, strand);

        let mut alignment = AlignmentRecord {
            name,
            reference_start,
            reference_end,
            sequence,
            quality_scores,
            mapping_quality,
            flags: u16::from(flags),
            strand,
            cigar: cigar_ops,
            variants: Vec::new(),
            methylation,
        };

        // Extract variants from CIGAR
        alignment.extract_variants();

        Ok(Some(alignment))
    }

    /// Clear all loaded data to free memory
    pub fn clear_all(&mut self)
    {
        self.loaded_tracks.clear();
    }

    /// Clear loaded data for specific chromosome
    pub fn clear_chromosome(&mut self, chromosome: &str)
    {
        self.loaded_tracks.remove(chromosome);
    }
}

/// Greedy row assignment - assign each read to first non-overlapping row
pub fn stack_alignments(
    records: &[AlignmentRecord],
    viewport_start: usize,
    viewport_end: usize,
    max_rows: usize,
) -> Vec<Vec<usize>>
{
    let mut rows: Vec<Vec<usize>> = Vec::new();

    for (idx, record) in records.iter().enumerate()
    {
        // Skip reads outside viewport
        if record.reference_end < viewport_start || record.reference_start > viewport_end
        {
            continue;
        }

        let mut placed = false;

        // Try to place in existing row
        for row in &mut rows
        {
            let can_place = row.iter().all(|&other_idx| {
                let other = &records[other_idx];
                record.reference_start >= other.reference_end
                    || record.reference_end <= other.reference_start
            });

            if can_place
            {
                row.push(idx);
                placed = true;
                break;
            }
        }

        // Create new row if needed and under limit
        if !placed && rows.len() < max_rows
        {
            rows.push(vec![idx]);
        }
    }

    rows
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn test_cigar_reference_length()
    {
        let cigar = vec![
            CigarOp::Match(10),
            CigarOp::Deletion(5),
            CigarOp::Insertion(3),
            CigarOp::Match(10),
        ];
        assert_eq!(AlignmentRecord::compute_reference_length(&cigar), 25);
    }

    #[test]
    fn test_stacking()
    {
        let records = vec![
            AlignmentRecord {
                name: "read1".to_string(),
                reference_start: 100,
                reference_end: 200,
                sequence: vec![],
                quality_scores: vec![],
                mapping_quality: 60,
                flags: 0,
                strand: Strand::Forward,
                cigar: vec![],
                variants: vec![],
                methylation: vec![],
            },
            AlignmentRecord {
                name: "read2".to_string(),
                reference_start: 150,
                reference_end: 250,
                sequence: vec![],
                quality_scores: vec![],
                mapping_quality: 60,
                flags: 0,
                strand: Strand::Forward,
                cigar: vec![],
                variants: vec![],
                methylation: vec![],
            },
            AlignmentRecord {
                name: "read3".to_string(),
                reference_start: 300,
                reference_end: 400,
                sequence: vec![],
                quality_scores: vec![],
                mapping_quality: 60,
                flags: 0,
                strand: Strand::Forward,
                cigar: vec![],
                variants: vec![],
                methylation: vec![],
            },
        ];

        let rows = stack_alignments(&records, 0, 500, 50);
        assert_eq!(rows.len(), 2); // read1 and read2 overlap, read3 doesn't
        assert!(rows[0].contains(&0) || rows[1].contains(&0));
        assert!(rows[0].contains(&1) || rows[1].contains(&1));
        assert!(rows[0].contains(&2) || rows[1].contains(&2));
    }
}
