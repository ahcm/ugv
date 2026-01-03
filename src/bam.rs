use crate::gff::Strand;
use anyhow::Result;
use std::collections::HashMap;
use std::io::{BufRead, Cursor};

use noodles::bam;
use noodles::sam;

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
        if self.reference_length == 0
        {
            return;
        }

        let num_bins = (self.reference_length + bin_size - 1) / bin_size;
        let mut coverage = vec![CoveragePoint::default(); num_bins];

        for record in &self.records
        {
            let start_bin = record.reference_start / bin_size;
            let end_bin = record.reference_end / bin_size;

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

                coverage[bin_idx].position = bin_idx * bin_size;
            }
        }

        self.coverage = coverage;
    }
}

/// Container for all tracks
#[derive(Clone)]
pub struct AlignmentData
{
    pub tracks: HashMap<String, AlignmentTrack>,
    pub header: sam::Header,
}

impl AlignmentData
{
    /// Load BAM from file (native only)
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(bam_path: &str) -> Result<Self>
    {
        use std::fs::File;

        let mut reader = bam::io::Reader::new(File::open(bam_path)?);
        let header = reader.read_header()?;

        Self::read_records(reader, header)
    }

    /// Load BAM from bytes (for WASM)
    pub fn from_bytes(data: Vec<u8>) -> Result<Self>
    {
        let cursor = Cursor::new(data);
        let mut reader = bam::io::Reader::new(cursor);
        let header = reader.read_header()?;

        Self::read_records(reader, header)
    }

    /// Read all records from BAM reader
    fn read_records<R: BufRead>(mut reader: bam::io::Reader<R>, header: sam::Header)
    -> Result<Self>
    {
        let mut tracks: HashMap<String, AlignmentTrack> = HashMap::new();

        // Initialize tracks from header
        for (name, reference_sequence) in header.reference_sequences().iter()
        {
            let name_str = name.to_string();
            let length: usize = reference_sequence.length().into();
            tracks.insert(name_str.clone(), AlignmentTrack::new(name_str, length));
        }

        // Read all records using the iterator
        for result in reader.records()
        {
            let record = result?;

            // Skip unmapped reads
            let flags = record.flags();
            if flags.is_unmapped()
            {
                continue;
            }

            // Get reference sequence ID
            let ref_seq_id = match record.reference_sequence_id()
            {
                Some(Ok(id)) => id,
                _ => continue,
            };

            // Get reference sequence name from header
            let ref_name = header
                .reference_sequences()
                .get_index(ref_seq_id)
                .map(|(name, _)| name.to_string());

            if ref_name.is_none()
            {
                continue;
            }
            let ref_name = ref_name.unwrap();

            // Get alignment start (convert from 1-based to 0-based)
            let reference_start = match record.alignment_start()
            {
                Some(Ok(pos)) => usize::from(pos) - 1,
                _ => continue,
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

            // Add to appropriate track
            if let Some(track) = tracks.get_mut(&ref_name)
            {
                track.records.push(alignment);
            }
        }

        // Compute coverage for all tracks
        for track in tracks.values_mut()
        {
            track.compute_coverage(100); // 100bp bins
        }

        Ok(AlignmentData { tracks, header })
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
