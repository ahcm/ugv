use crate::gff::Strand;
use anyhow::Result;
use std::collections::HashMap;

use noodles::bam;
use noodles::sam;

#[cfg(target_arch = "wasm32")]
use std::io::Cursor;

// Memory management constants
const MAX_LOADED_REGION_SIZE: usize = 10_000_000; // 10MB max region
const MAX_LOADED_ALIGNMENTS: usize = 100_000; // Max 100k reads
const VIEWPORT_BUFFER: usize = 2_000_000; // ±2Mb buffer around viewport

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

/// Single alignment record with pre-computed fields
#[derive(Debug, Clone)]
pub struct AlignmentRecord
{
    pub name: String,
    pub reference_start: usize,  // 0-based position
    pub reference_end: usize,    // Computed from alignment length
    pub sequence: Vec<u8>,       // Read sequence
    pub quality_scores: Vec<u8>, // Phred quality scores
    pub mapping_quality: u8,     // MAPQ
    pub flags: u16,              // SAM flags
    pub strand: Strand,          // Forward or Reverse
    pub cigar: Vec<CigarOp>,     // Parsed CIGAR operations
    pub variants: Vec<Variant>,  // Pre-computed from CIGAR + MD tag
}

impl AlignmentRecord
{
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
        let region_length = region_end - region_start;
        let num_bins = (region_length + bin_size - 1) / bin_size;
        let mut coverage = vec![CoveragePoint::default(); num_bins];

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

                    coverage[bin_idx].position = region_start + bin_idx * bin_size;
                }
            }
        }

        self.coverage = coverage;
    }

    /// Clear loaded data to free memory
    pub fn clear(&mut self)
    {
        self.records.clear();
        self.coverage.clear();
        self.loaded_region = None;
    }
}

/// Container for BAM data with lazy loading support
#[derive(Clone)]
pub struct AlignmentData
{
    #[cfg(not(target_arch = "wasm32"))]
    pub bam_path: String,

    #[cfg(target_arch = "wasm32")]
    pub bam_bytes: Vec<u8>,

    pub header: sam::Header,
    pub reference_lengths: HashMap<String, usize>,

    // Currently loaded tracks (only for displayed region)
    pub loaded_tracks: HashMap<String, AlignmentTrack>,
}

impl AlignmentData
{
    /// Create from file path (native only) - does NOT load records yet
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(bam_path: &str) -> Result<Self>
    {
        use std::fs::File;

        let mut reader = bam::io::Reader::new(File::open(bam_path)?);
        let header = reader.read_header()?;

        // Extract reference lengths from header
        let mut reference_lengths = HashMap::new();
        for (name, reference_sequence) in header.reference_sequences().iter()
        {
            let name_str = name.to_string();
            let length: usize = reference_sequence.length().into();
            reference_lengths.insert(name_str, length);
        }

        Ok(AlignmentData {
            bam_path: bam_path.to_string(),
            header,
            reference_lengths,
            loaded_tracks: HashMap::new(),
        })
    }

    /// Create from bytes (WASM) - does NOT load records yet
    #[cfg(target_arch = "wasm32")]
    pub fn from_bytes(data: Vec<u8>) -> Result<Self>
    {
        let cursor = Cursor::new(&data);
        let mut reader = bam::io::Reader::new(cursor);
        let header = reader.read_header()?;

        // Extract reference lengths from header
        let mut reference_lengths = HashMap::new();
        for (name, reference_sequence) in header.reference_sequences().iter()
        {
            let name_str = name.to_string();
            let length: usize = reference_sequence.length().into();
            reference_lengths.insert(name_str, length);
        }

        Ok(AlignmentData {
            bam_bytes: data,
            header,
            reference_lengths,
            loaded_tracks: HashMap::new(),
        })
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
        // Calculate region with buffer
        let buffer_start = viewport_start.saturating_sub(VIEWPORT_BUFFER);
        let buffer_end = (viewport_end + VIEWPORT_BUFFER).min(
            *self.reference_lengths.get(chromosome).unwrap_or(&usize::MAX)
        );

        // Limit region size
        let region_size = buffer_end - buffer_start;
        if region_size > MAX_LOADED_REGION_SIZE
        {
            // Region too large, load only viewport without buffer
            return self.query_region_exact(chromosome, viewport_start, viewport_end);
        }

        // Check if we need to load (avoid borrow checker issues)
        let need_load = if let Some(track) = self.loaded_tracks.get(chromosome)
        {
            if let Some((loaded_start, loaded_end)) = track.loaded_region
            {
                // Check if requested region is outside loaded region
                buffer_start < loaded_start || buffer_end > loaded_end
            }
            else
            {
                true
            }
        }
        else
        {
            true
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
    fn load_region(
        &mut self,
        chromosome: &str,
        start: usize,
        end: usize,
    ) -> Result<()>
    {
        use noodles::core::Region;

        let reference_length = *self.reference_lengths.get(chromosome).unwrap_or(&0);
        let mut track = AlignmentTrack::new(chromosome.to_string(), reference_length);
        track.loaded_region = Some((start, end));

        // Try to use indexed reader
        let indexed_result = bam::io::indexed_reader::Builder::default()
            .build_from_path(&self.bam_path);

        if let Ok(mut reader) = indexed_result
        {
            // Use indexed query for efficient region loading
            let region = format!("{}:{}-{}", chromosome, start + 1, end); // 1-based for noodles
            if let Ok(region) = region.parse::<Region>()
            {
                let header = reader.read_header()?;
                let query = reader.query(&header, &region)?;

                let mut count = 0;
                for result in query
                {
                    if count >= MAX_LOADED_ALIGNMENTS
                    {
                        eprintln!("Warning: Hit alignment limit of {} for region", MAX_LOADED_ALIGNMENTS);
                        break;
                    }

                    let record = result?;
                    if let Some(alignment) = Self::parse_record(&record, &header)?
                    {
                        track.records.push(alignment);
                        count += 1;
                    }
                }
            }
        }
        else
        {
            // Fallback: no index available, scan entire file (slow)
            eprintln!("Warning: BAM index not found for {}, scanning entire file (slow)", self.bam_path);

            use std::fs::File;
            let mut reader = bam::io::Reader::new(File::open(&self.bam_path)?);
            let header = reader.read_header()?;

            let mut count = 0;
            for result in reader.records()
            {
                if count >= MAX_LOADED_ALIGNMENTS
                {
                    break;
                }

                let record = result?;

                // Check if record is in our region
                if let Some(Ok(ref_seq_id)) = record.reference_sequence_id()
                {
                    if let Some((ref_name, _)) = header.reference_sequences().get_index(ref_seq_id)
                    {
                        if ref_name.to_string() == chromosome
                        {
                            if let Some(Ok(pos)) = record.alignment_start()
                            {
                                let pos_0based = usize::from(pos) - 1;
                                if pos_0based >= start && pos_0based <= end
                                {
                                    if let Some(alignment) = Self::parse_record(&record, &header)?
                                    {
                                        track.records.push(alignment);
                                        count += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Compute coverage for this region
        track.compute_coverage(100); // 100bp bins

        // Store in loaded_tracks, replacing any previous data for this chromosome
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

        let reference_length = *self.reference_lengths.get(chromosome).unwrap_or(&0);
        let mut track = AlignmentTrack::new(chromosome.to_string(), reference_length);
        track.loaded_region = Some((0, reference_length));

        let mut count = 0;
        for result in reader.records()
        {
            if count >= MAX_LOADED_ALIGNMENTS
            {
                eprintln!("Warning: Hit alignment limit of {} for chromosome", MAX_LOADED_ALIGNMENTS);
                break;
            }

            let record = result?;

            // Check if this record belongs to our chromosome
            if let Some(Ok(ref_seq_id)) = record.reference_sequence_id()
            {
                if let Some((ref_name, _)) = header.reference_sequences().get_index(ref_seq_id)
                {
                    if ref_name.to_string() == chromosome
                    {
                        if let Some(alignment) = Self::parse_record(&record, &header)?
                        {
                            track.records.push(alignment);
                            count += 1;
                        }
                    }
                }
            }
        }

        track.compute_coverage(100);
        self.loaded_tracks.insert(chromosome.to_string(), track);

        Ok(())
    }

    /// Parse a single BAM record into AlignmentRecord
    fn parse_record(
        record: &bam::Record,
        header: &sam::Header,
    ) -> Result<Option<AlignmentRecord>>
    {
        // Skip unmapped reads
        let flags = record.flags();
        if flags.is_unmapped()
        {
            return Ok(None);
        }

        // Get reference sequence ID
        let ref_seq_id = match record.reference_sequence_id()
        {
            Some(Ok(id)) => id,
            _ => return Ok(None),
        };

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
        let sequence: Vec<u8> = (0..record.sequence().len())
            .filter_map(|i| record.sequence().get(i))
            .collect();

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
            },
        ];

        let rows = stack_alignments(&records, 0, 500, 50);
        assert_eq!(rows.len(), 2); // read1 and read2 overlap, read3 doesn't
        assert!(rows[0].contains(&0) || rows[1].contains(&0));
        assert!(rows[0].contains(&1) || rows[1].contains(&1));
        assert!(rows[0].contains(&2) || rows[1].contains(&2));
    }
}
