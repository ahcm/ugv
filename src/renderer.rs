use crate::bam::{AlignmentRecord, AlignmentTrack, VariantType};
use crate::fasta::Chromosome;
use crate::gff::Feature;
use crate::interval_tree::IntervalTree;
use crate::translation;
use crate::tsv::TsvChromosomeTrack;
use crate::viewport::Viewport;
use egui::{Color32, FontId, Painter, Pos2, Rect, Stroke, Vec2};

const RULER_HEIGHT: f32 = 30.0;
const SEQUENCE_HEIGHT: f32 = 40.0;
const GC_CONTENT_HEIGHT: f32 = 60.0;
const FEATURE_TRACK_HEIGHT: f32 = 30.0;
const AMINO_ACID_TRACK_HEIGHT: f32 = 20.0;
const TRACK_SPACING: f32 = 10.0;

// BAM track constants
const COVERAGE_TRACK_HEIGHT: f32 = 80.0;
const ALIGNMENT_ROW_HEIGHT: f32 = 12.0;
const MAX_ALIGNMENT_ROWS: usize = 50;
const VARIANT_TRACK_HEIGHT: f32 = 30.0;

// TSV track constants
const TSV_TRACK_HEIGHT: f32 = 100.0;

// Helper function to draw an empty track placeholder
pub fn draw_empty_track(
    painter: &Painter,
    rect: Rect,
    y_offset: f32,
    height: f32,
    label: &str,
) -> f32
{
    let track_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));

    // Background
    painter.rect_filled(track_rect, 0.0, Color32::from_gray(250));

    // Border
    painter.rect_stroke(track_rect, 0.0, Stroke::new(1.0, Color32::from_gray(220)));

    // Draw label
    painter.text(
        Pos2::new(rect.left() + 5.0, y_offset + 5.0),
        egui::Align2::LEFT_TOP,
        label,
        FontId::proportional(12.0),
        Color32::from_gray(120),
    );

    // Draw "No data" message in center
    painter.text(
        Pos2::new(rect.center().x, y_offset + height / 2.0),
        egui::Align2::CENTER_CENTER,
        "No data in viewport",
        FontId::proportional(11.0),
        Color32::from_gray(150),
    );

    track_rect.bottom()
}

// Make these functions public so they can be called from main.rs
pub fn draw_ruler(
    painter: &Painter,
    rect: Rect,
    viewport: &Viewport,
    y_offset: f32,
    height: f32,
    chromosome_name: &str,
    selection: Option<(usize, usize)>,
) -> f32
{
    let ruler_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));

    // Background
    painter.rect_filled(ruler_rect, 0.0, Color32::from_gray(240));

    if let Some((selection_start, selection_end)) = selection
    {
        let start = selection_start.min(selection_end);
        let end = selection_start.max(selection_end).saturating_add(1);
        let clipped_start = start.max(viewport.start).min(viewport.end);
        let clipped_end = end.max(viewport.start).min(viewport.end);

        if clipped_end > clipped_start
        {
            let x_start = viewport.position_to_screen(clipped_start, rect.width()) + rect.left();
            let x_end = viewport.position_to_screen(clipped_end, rect.width()) + rect.left();
            let selection_rect = Rect::from_min_max(
                Pos2::new(x_start, ruler_rect.top()),
                Pos2::new(x_end, ruler_rect.bottom()),
            );
            painter.rect_filled(
                selection_rect,
                0.0,
                Color32::from_rgba_unmultiplied(120, 170, 220, 120),
            );
            painter.rect_stroke(
                selection_rect,
                0.0,
                Stroke::new(1.0, Color32::from_rgb(80, 120, 170)),
            );
        }
    }

    // Draw chromosome name on the left
    painter.text(
        Pos2::new(rect.left() + 10.0, ruler_rect.top() + height / 2.0),
        egui::Align2::LEFT_CENTER,
        chromosome_name,
        egui::FontId::proportional(14.0),
        Color32::from_rgb(60, 60, 100),
    );

    // Draw tick marks
    let width = viewport.width();
    let tick_interval = calculate_tick_interval(width);

    let first_tick = ((viewport.start / tick_interval) + 1) * tick_interval;

    for tick in (first_tick..viewport.end).step_by(tick_interval)
    {
        let x = viewport.position_to_screen(tick, rect.width()) + rect.left();

        // Tick line
        painter.line_segment(
            [
                Pos2::new(x, ruler_rect.top()),
                Pos2::new(x, ruler_rect.top() + 10.0),
            ],
            Stroke::new(1.0, Color32::DARK_GRAY),
        );

        // Label
        let label = format_position(tick, tick_interval);
        painter.text(
            Pos2::new(x, ruler_rect.top() + 15.0),
            egui::Align2::CENTER_TOP,
            label,
            egui::FontId::monospace(10.0),
            Color32::BLACK,
        );
    }

    ruler_rect.bottom()
}

pub fn draw_gc_content(
    painter: &Painter,
    rect: Rect,
    chromosome: &Chromosome,
    viewport: &Viewport,
    y_offset: f32,
    height: f32,
) -> f32
{
    let gc_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));

    // Background
    painter.rect_filled(gc_rect, 0.0, Color32::from_gray(250));

    // Draw 50% reference line
    let mid_y = gc_rect.center().y;
    painter.line_segment(
        [
            Pos2::new(gc_rect.left(), mid_y),
            Pos2::new(gc_rect.right(), mid_y),
        ],
        Stroke::new(1.0, Color32::from_gray(200)),
    );

    // Calculate GC content in windows
    let window_size = (viewport.width() / rect.width() as usize).max(100);
    let mut points = Vec::new();

    for i in (viewport.start..viewport.end).step_by(window_size.max(1))
    {
        let gc = chromosome.get_gc_content(i, (i + window_size).min(viewport.end));
        let x = viewport.position_to_screen(i, rect.width()) + rect.left();
        let y = gc_rect.bottom() - (gc * height);
        points.push(Pos2::new(x, y));
    }

    // Draw GC content line
    if points.len() > 1
    {
        for window in points.windows(2)
        {
            painter.line_segment(
                [window[0], window[1]],
                Stroke::new(2.0, Color32::from_rgb(70, 130, 180)),
            );
        }
    }

    // Label
    painter.text(
        Pos2::new(gc_rect.left() + 5.0, gc_rect.top() + 5.0),
        egui::Align2::LEFT_TOP,
        "GC%",
        egui::FontId::monospace(10.0),
        Color32::DARK_GRAY,
    );

    gc_rect.bottom()
}

pub fn draw_sequence(
    painter: &Painter,
    rect: Rect,
    chromosome: &Chromosome,
    viewport: &Viewport,
    y_offset: f32,
    height: f32,
) -> f32
{
    let seq_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));

    painter.rect_filled(seq_rect, 0.0, Color32::from_gray(245));

    let pixels_per_base = rect.width() / viewport.width() as f32;

    if pixels_per_base >= 8.0
    {
        // Draw individual bases
        for pos in viewport.start..viewport.end.min(chromosome.length)
        {
            if let Some(base) = chromosome.get_base_at(pos)
            {
                let x = viewport.position_to_screen(pos, rect.width()) + rect.left();

                let color = match base
                {
                    'A' => Color32::from_rgb(0, 200, 0),
                    'T' => Color32::from_rgb(200, 0, 0),
                    'G' => Color32::from_rgb(0, 0, 200),
                    'C' => Color32::from_rgb(200, 200, 0),
                    _ => Color32::GRAY,
                };

                painter.text(
                    Pos2::new(x, seq_rect.center().y),
                    egui::Align2::CENTER_CENTER,
                    base.to_string(),
                    egui::FontId::monospace(10.0),
                    color,
                );
            }
        }
    }
    else
    {
        // Draw colored bars representing base composition
        for pos in viewport.start..viewport.end.min(chromosome.length)
        {
            if let Some(base) = chromosome.get_base_at(pos)
            {
                let x = viewport.position_to_screen(pos, rect.width()) + rect.left();

                let color = match base
                {
                    'A' => Color32::from_rgb(0, 200, 0),
                    'T' => Color32::from_rgb(200, 0, 0),
                    'G' => Color32::from_rgb(0, 0, 200),
                    'C' => Color32::from_rgb(200, 200, 0),
                    _ => Color32::GRAY,
                };

                painter.rect_filled(
                    Rect::from_min_size(
                        Pos2::new(x, seq_rect.top()),
                        Vec2::new(pixels_per_base.max(1.0), SEQUENCE_HEIGHT),
                    ),
                    0.0,
                    color,
                );
            }
        }
    }

    seq_rect.bottom()
}

pub fn draw_amino_acid_frames(
    painter: &Painter,
    rect: Rect,
    chromosome: &Chromosome,
    viewport: &Viewport,
    y_offset: f32,
    height: f32,
) -> f32
{
    // Get the visible sequence
    let start = viewport.start;
    let end = viewport.end.min(chromosome.length);

    if start >= end
    {
        return y_offset;
    }

    let sequence = &chromosome.sequence[start..end];
    let sequence_length = sequence.len();

    // Translate all 6 frames
    let (forward_frames, reverse_frames) = translation::translate_six_frames(sequence);

    // Divide height equally among 6 frames
    let frame_height = height / 6.0;

    let frame_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));

    // Background
    painter.rect_filled(frame_rect, 0.0, Color32::from_gray(250));

    let mut current_y = y_offset;

    // Draw forward frames (top 3)
    for (frame_idx, amino_acids) in forward_frames.iter().enumerate()
    {
        draw_single_frame(
            painter,
            rect,
            amino_acids,
            viewport,
            start,
            sequence_length,
            frame_idx,
            current_y,
            frame_height,
            true, // forward strand
        );
        current_y += frame_height;
    }

    // Draw reverse frames (bottom 3)
    for (frame_idx, amino_acids) in reverse_frames.iter().enumerate()
    {
        draw_single_frame(
            painter,
            rect,
            amino_acids,
            viewport,
            start,
            sequence_length,
            frame_idx,
            current_y,
            frame_height,
            false, // reverse strand
        );
        current_y += frame_height;
    }

    frame_rect.bottom()
}

fn draw_single_frame(
    painter: &Painter,
    rect: Rect,
    amino_acids: &[char],
    viewport: &Viewport,
    sequence_start: usize,
    sequence_length: usize,
    frame: usize,
    y_offset: f32,
    height: f32,
    is_forward: bool,
)
{
    let track_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));

    // Draw background
    let bg_color = if is_forward
    {
        Color32::from_rgb(240, 248, 255) // Light blue for forward
    }
    else
    {
        Color32::from_rgb(255, 248, 240) // Light orange for reverse
    };
    painter.rect_filled(track_rect, 0.0, bg_color);

    // Draw border
    painter.rect_stroke(track_rect, 0.0, Stroke::new(1.0, Color32::from_gray(200)));

    // Draw frame label
    let strand_symbol = if is_forward { "+" } else { "-" };
    let frame_label = format!("{}{}", strand_symbol, frame + 1);
    painter.text(
        Pos2::new(rect.left() + 5.0, y_offset + AMINO_ACID_TRACK_HEIGHT / 2.0),
        egui::Align2::LEFT_CENTER,
        frame_label,
        FontId::monospace(10.0),
        Color32::from_gray(100),
    );

    // Draw amino acids
    let pixels_per_base = rect.width() / viewport.width() as f32;
    let char_width = pixels_per_base * 3.0; // Each amino acid represents 3 bases (codon)

    // Only draw if characters are visible
    if char_width > 6.0
    {
        for (aa_idx, &aa) in amino_acids.iter().enumerate()
        {
            let codon_pos = if is_forward
            {
                sequence_start + frame + (aa_idx * 3)
            }
            else
            {
                // For reverse strand, calculate the position
                sequence_start + sequence_length.saturating_sub(frame + (aa_idx * 3) + 3)
            };

            if codon_pos >= viewport.start && codon_pos < viewport.end
            {
                let x = viewport.position_to_screen(codon_pos, rect.width()) + rect.left();

                // Color based on amino acid type
                let color = get_amino_acid_color(aa);

                // Draw background for the amino acid
                if char_width > 10.0
                {
                    painter.rect_filled(
                        Rect::from_min_size(
                            Pos2::new(x, track_rect.top()),
                            Vec2::new(char_width.min(rect.right() - x), AMINO_ACID_TRACK_HEIGHT),
                        ),
                        0.0,
                        color.linear_multiply(0.3),
                    );
                }

                // Draw the amino acid letter if there's enough space
                if char_width > 12.0
                {
                    painter.text(
                        Pos2::new(x + char_width / 2.0, y_offset + AMINO_ACID_TRACK_HEIGHT / 2.0),
                        egui::Align2::CENTER_CENTER,
                        aa.to_string(),
                        FontId::monospace(10.0),
                        color,
                    );
                }
            }
        }
    }
}

fn get_amino_acid_color(aa: char) -> Color32
{
    match aa
    {
        // Hydrophobic (aliphatic)
        'A' | 'V' | 'I' | 'L' | 'M' => Color32::from_rgb(0, 150, 0),
        // Aromatic
        'F' | 'W' | 'Y' => Color32::from_rgb(100, 50, 150),
        // Polar
        'S' | 'T' | 'N' | 'Q' => Color32::from_rgb(0, 100, 150),
        // Positively charged
        'K' | 'R' | 'H' => Color32::from_rgb(0, 0, 200),
        // Negatively charged
        'D' | 'E' => Color32::from_rgb(200, 0, 0),
        // Special
        'C' => Color32::from_rgb(200, 200, 0), // Cysteine (yellow)
        'G' => Color32::from_rgb(150, 150, 150), // Glycine (gray)
        'P' => Color32::from_rgb(200, 100, 0), // Proline (orange)
        // Stop codon
        '*' => Color32::from_rgb(200, 0, 0),
        // Unknown
        _ => Color32::GRAY,
    }
}

pub fn draw_features(
    painter: &Painter,
    rect: Rect,
    features: &[Feature],
    viewport: &Viewport,
    interval_tree: &IntervalTree,
    chr_name: &str,
    y_offset: f32,
    height: f32,
) -> f32
{
    // Query features in viewport
    let feature_indices = interval_tree.query(chr_name, viewport.start, viewport.end);

    if feature_indices.is_empty()
    {
        return draw_empty_track(painter, rect, y_offset, height, "Features");
    }

    // Group features by type for different tracks
    let mut tracks: Vec<Vec<usize>> = Vec::new();

    for &idx in &feature_indices
    {
        let feature = &features[idx];

        // Find a track where this feature fits
        let mut placed = false;
        for track in &mut tracks
        {
            let can_place = track.iter().all(|&other_idx| {
                let other = &features[other_idx];
                feature.start >= other.end || feature.end <= other.start
            });

            if can_place
            {
                track.push(idx);
                placed = true;
                break;
            }
        }

        if !placed
        {
            tracks.push(vec![idx]);
        }
    }

    // Calculate row height based on total height and number of rows
    let num_rows = tracks.len().max(1);
    let row_height = if num_rows > 1
    {
        (height - (TRACK_SPACING * (num_rows - 1) as f32)) / num_rows as f32
    }
    else
    {
        height
    };

    // Draw each track
    for (track_idx, track) in tracks.iter().enumerate()
    {
        let track_y = y_offset + (track_idx as f32 * (row_height + TRACK_SPACING));

        for &idx in track
        {
            let feature = &features[idx];
            draw_feature(painter, rect, feature, viewport, track_y, row_height);
        }
    }

    y_offset + height
}

fn draw_feature(
    painter: &Painter,
    rect: Rect,
    feature: &Feature,
    viewport: &Viewport,
    y_offset: f32,
    height: f32,
)
{
    let start_x = viewport.position_to_screen(feature.start, rect.width()) + rect.left();
    let end_x = viewport.position_to_screen(feature.end, rect.width()) + rect.left();

    let feature_rect = Rect::from_min_size(
        Pos2::new(start_x, y_offset),
        Vec2::new((end_x - start_x).max(2.0), height),
    );

    // Draw feature box
    painter.rect_filled(feature_rect, 2.0, feature.color_by_type());
    painter.rect_stroke(feature_rect, 2.0, Stroke::new(1.0, Color32::BLACK));

    // Draw strand indicator
    match feature.strand
    {
        crate::gff::Strand::Forward =>
        {
            let arrow_x = feature_rect.right() - 5.0;
            painter.line_segment(
                [
                    Pos2::new(arrow_x - 5.0, feature_rect.center().y - 3.0),
                    Pos2::new(arrow_x, feature_rect.center().y),
                ],
                Stroke::new(2.0, Color32::BLACK),
            );
            painter.line_segment(
                [
                    Pos2::new(arrow_x - 5.0, feature_rect.center().y + 3.0),
                    Pos2::new(arrow_x, feature_rect.center().y),
                ],
                Stroke::new(2.0, Color32::BLACK),
            );
        }
        crate::gff::Strand::Reverse =>
        {
            let arrow_x = feature_rect.left() + 5.0;
            painter.line_segment(
                [
                    Pos2::new(arrow_x + 5.0, feature_rect.center().y - 3.0),
                    Pos2::new(arrow_x, feature_rect.center().y),
                ],
                Stroke::new(2.0, Color32::BLACK),
            );
            painter.line_segment(
                [
                    Pos2::new(arrow_x + 5.0, feature_rect.center().y + 3.0),
                    Pos2::new(arrow_x, feature_rect.center().y),
                ],
                Stroke::new(2.0, Color32::BLACK),
            );
        }
        _ =>
        {}
    }

    // Draw label if there's enough space
    if feature_rect.width() > 50.0
    {
        painter.text(
            feature_rect.center(),
            egui::Align2::CENTER_CENTER,
            feature.name(),
            egui::FontId::monospace(9.0),
            Color32::WHITE,
        );
    }
}

fn calculate_tick_interval(width: usize) -> usize
{
    let base_intervals = [
        1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10_000, 20_000, 50_000, 100_000,
        200_000, 500_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000,
    ];

    for &interval in &base_intervals
    {
        if width / interval < 10
        {
            return interval;
        }
    }

    10_000_000
}

fn format_position(pos: usize, tick_interval: usize) -> String
{
    let format_with_unit = |value: f64, unit: &str, interval_value: f64| {
        let decimals = if interval_value >= 1.0
        {
            0
        }
        else if interval_value >= 0.1
        {
            1
        }
        else if interval_value >= 0.01
        {
            2
        }
        else
        {
            3
        };
        if decimals == 0
        {
            format!("{:.0}{}", value, unit)
        }
        else
        {
            format!("{:.*}{}", decimals, value, unit)
        }
    };

    if pos >= 1_000_000
    {
        let scale = 1_000_000.0;
        let value = pos as f64 / scale;
        let interval_value = tick_interval as f64 / scale;
        format_with_unit(value, "Mb", interval_value)
    }
    else if pos >= 1_000
    {
        let scale = 1_000.0;
        let value = pos as f64 / scale;
        let interval_value = tick_interval as f64 / scale;
        format_with_unit(value, "kb", interval_value)
    }
    else
    {
        pos.to_string()
    }
}

// ============================================================================
// BAM TRACK RENDERING
// ============================================================================

/// Draw coverage histogram track
pub fn draw_coverage_track(
    painter: &Painter,
    rect: Rect,
    track: &AlignmentTrack,
    viewport: &Viewport,
    y_offset: f32,
    height: f32,
) -> f32
{
    let track_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));

    // Background
    painter.rect_filled(track_rect, 0.0, Color32::from_gray(250));

    // Draw label
    painter.text(
        Pos2::new(rect.left() + 5.0, y_offset + 5.0),
        egui::Align2::LEFT_TOP,
        "Coverage",
        FontId::proportional(12.0),
        Color32::BLACK,
    );

    // Get coverage data for viewport
    let bin_size = 100; // matches the bin size used in coverage calculation

    // Check if we have a loaded region and coverage data
    if track.coverage.is_empty()
    {
        painter.text(
            Pos2::new(rect.center().x, y_offset + height / 2.0),
            egui::Align2::CENTER_CENTER,
            "No coverage data",
            FontId::proportional(11.0),
            Color32::from_gray(150),
        );
        return track_rect.bottom();
    }

    // Calculate bin indices relative to loaded region
    let (region_start, region_end) = track.loaded_region.unwrap_or((0, track.reference_length));

    // Only show coverage if viewport overlaps with loaded region
    if viewport.end <= region_start || viewport.start >= region_end
    {
        painter.text(
            Pos2::new(rect.center().x, y_offset + height / 2.0),
            egui::Align2::CENTER_CENTER,
            "Coverage not loaded for this region",
            FontId::proportional(11.0),
            Color32::from_gray(150),
        );
        return track_rect.bottom();
    }

    // Calculate bin indices relative to the loaded region start
    let viewport_start_in_region = viewport.start.saturating_sub(region_start);
    let viewport_end_in_region = viewport.end.saturating_sub(region_start);

    let start_bin = viewport_start_in_region / bin_size;
    let end_bin = ((viewport_end_in_region + bin_size - 1) / bin_size).min(track.coverage.len());

    if start_bin >= track.coverage.len()
    {
        painter.text(
            Pos2::new(rect.center().x, y_offset + height / 2.0),
            egui::Align2::CENTER_CENTER,
            "No coverage in this region",
            FontId::proportional(11.0),
            Color32::from_gray(150),
        );
        return track_rect.bottom();
    }

    let coverage_slice = &track.coverage[start_bin..end_bin];

    if coverage_slice.is_empty()
    {
        painter.text(
            Pos2::new(rect.center().x, y_offset + height / 2.0),
            egui::Align2::CENTER_CENTER,
            "No coverage in this region",
            FontId::proportional(11.0),
            Color32::from_gray(150),
        );
        return track_rect.bottom();
    }

    // Find max depth for scaling
    let max_depth = coverage_slice
        .iter()
        .map(|p| p.depth)
        .max()
        .unwrap_or(1)
        .max(1);

    // Draw histogram bars
    for point in coverage_slice
    {
        if point.position < viewport.start || point.position >= viewport.end
        {
            continue;
        }

        let x = viewport.position_to_screen(point.position, rect.width()) + rect.left();
        let bar_width = (rect.width() / viewport.width() as f32 * bin_size as f32).max(1.0);
        let bar_height = (point.depth as f32 / max_depth as f32) * (height - 20.0);
        let y = track_rect.bottom() - bar_height - 5.0;

        let color = Color32::from_rgb(100, 150, 200);
        painter.rect_filled(
            Rect::from_min_size(Pos2::new(x, y), Vec2::new(bar_width, bar_height)),
            0.0,
            color,
        );
    }

    // Update label with max depth
    let label_text = format!("Coverage (max: {})", max_depth);
    painter.text(
        Pos2::new(rect.left() + 5.0, y_offset + 5.0),
        egui::Align2::LEFT_TOP,
        label_text,
        FontId::proportional(12.0),
        Color32::BLACK,
    );

    track_rect.bottom()
}

/// Draw alignment reads track
pub fn draw_alignments(
    painter: &Painter,
    rect: Rect,
    track: &AlignmentTrack,
    viewport: &Viewport,
    y_offset: f32,
    height: f32,
    max_reads_display: usize,
) -> f32
{
    // Filter reads in viewport
    let visible_reads: Vec<&AlignmentRecord> = track
        .records
        .iter()
        .filter(|r| r.reference_end > viewport.start && r.reference_start < viewport.end)
        .collect();

    // If too many reads, show message instead
    if visible_reads.len() > max_reads_display
    {
        painter.text(
            Pos2::new(rect.left() + 5.0, y_offset + 5.0),
            egui::Align2::LEFT_TOP,
            format!("Too many reads ({}) - zoom in or increase limit to view", visible_reads.len()),
            FontId::proportional(12.0),
            Color32::DARK_GRAY,
        );
        return y_offset + height;
    }

    if visible_reads.is_empty()
    {
        return y_offset;
    }

    // Stack alignments
    let rows = crate::bam::stack_alignments(
        &track.records,
        viewport.start,
        viewport.end,
        MAX_ALIGNMENT_ROWS,
    );

    // Calculate row height based on total height
    let label_height = 20.0;
    let available_height = height - label_height;
    let num_rows = rows.len().min(MAX_ALIGNMENT_ROWS);
    let row_spacing = 2.0;
    let row_height = if num_rows > 0
    {
        ((available_height - (num_rows as f32 - 1.0) * row_spacing) / num_rows as f32).max(1.0)
    }
    else
    {
        12.0
    };

    // Background
    let track_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));
    painter.rect_filled(track_rect, 0.0, Color32::from_gray(245));

    // Draw label
    painter.text(
        Pos2::new(rect.left() + 5.0, y_offset + 5.0),
        egui::Align2::LEFT_TOP,
        format!("Alignments ({} reads)", visible_reads.len()),
        FontId::proportional(12.0),
        Color32::BLACK,
    );

    // Draw reads
    for (row_idx, row) in rows.iter().enumerate()
    {
        let row_y = y_offset + label_height + (row_idx as f32 * (row_height + row_spacing));

        for &read_idx in row
        {
            if read_idx < track.records.len()
            {
                let read = &track.records[read_idx];
                draw_single_read(painter, &rect, read, viewport, row_y, row_height);
            }
        }
    }

    track_rect.bottom()
}

/// Draw a single alignment read
fn draw_single_read(
    painter: &Painter,
    rect: &Rect,
    read: &AlignmentRecord,
    viewport: &Viewport,
    y: f32,
    height: f32,
)
{
    let start_x = viewport.position_to_screen(read.reference_start, rect.width()) + rect.left();
    let end_x = viewport.position_to_screen(read.reference_end, rect.width()) + rect.left();
    let width = (end_x - start_x).max(2.0);

    // Color by strand
    let base_color = match read.strand
    {
        crate::gff::Strand::Forward => Color32::from_rgb(100, 150, 200),
        crate::gff::Strand::Reverse => Color32::from_rgb(200, 100, 100),
        _ => Color32::GRAY,
    };

    // Adjust alpha by mapping quality
    let alpha = ((read.mapping_quality as f32 / 60.0) * 200.0 + 55.0).min(255.0) as u8;
    let color =
        Color32::from_rgba_premultiplied(base_color.r(), base_color.g(), base_color.b(), alpha);

    let read_rect = Rect::from_min_size(Pos2::new(start_x, y), Vec2::new(width, height));

    painter.rect_filled(read_rect, 1.0, color);
    painter.rect_stroke(read_rect, 1.0, Stroke::new(0.5, Color32::from_gray(100)));
}

/// Draw variant summary track
pub fn draw_variant_summary(
    painter: &Painter,
    rect: Rect,
    track: &AlignmentTrack,
    viewport: &Viewport,
    y_offset: f32,
    height: f32,
) -> f32
{
    let track_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));

    // Background
    painter.rect_filled(track_rect, 0.0, Color32::from_gray(250));

    // Collect variant counts by position
    let mut variant_counts: std::collections::HashMap<usize, (u32, u32, u32)> =
        std::collections::HashMap::new();

    for record in &track.records
    {
        if record.reference_start > viewport.end || record.reference_end < viewport.start
        {
            continue;
        }

        for variant in &record.variants
        {
            if variant.position >= viewport.start && variant.position < viewport.end
            {
                let counts = variant_counts.entry(variant.position).or_insert((0, 0, 0));
                match variant.variant_type
                {
                    VariantType::SNP => counts.0 += 1,
                    VariantType::Insertion => counts.1 += 1,
                    VariantType::Deletion => counts.2 += 1,
                }
            }
        }
    }

    // Draw label
    painter.text(
        Pos2::new(rect.left() + 5.0, y_offset + 5.0),
        egui::Align2::LEFT_TOP,
        format!("Variants ({} positions)", variant_counts.len()),
        FontId::proportional(12.0),
        Color32::BLACK,
    );

    // Draw variant markers
    for (pos, (snps, ins, del)) in variant_counts
    {
        let x = viewport.position_to_screen(pos, rect.width()) + rect.left();
        let total = snps + ins + del;
        let bar_height = ((total as f32 / 10.0).min(1.0) * (height - 15.0)).max(2.0);
        let y_bottom = track_rect.bottom() - 5.0;

        // Draw stacked bars
        let mut current_y = y_bottom;

        if del > 0
        {
            let del_height = (del as f32 / total as f32) * bar_height;
            painter.rect_filled(
                Rect::from_min_size(
                    Pos2::new(x - 1.0, current_y - del_height),
                    Vec2::new(2.0, del_height),
                ),
                0.0,
                Color32::BLACK,
            );
            current_y -= del_height;
        }

        if ins > 0
        {
            let ins_height = (ins as f32 / total as f32) * bar_height;
            painter.rect_filled(
                Rect::from_min_size(
                    Pos2::new(x - 1.0, current_y - ins_height),
                    Vec2::new(2.0, ins_height),
                ),
                0.0,
                Color32::from_rgb(150, 50, 200),
            );
            current_y -= ins_height;
        }

        if snps > 0
        {
            let snp_height = (snps as f32 / total as f32) * bar_height;
            painter.rect_filled(
                Rect::from_min_size(
                    Pos2::new(x - 1.0, current_y - snp_height),
                    Vec2::new(2.0, snp_height),
                ),
                0.0,
                Color32::from_rgb(200, 50, 50),
            );
        }
    }

    track_rect.bottom()
}

/// Draw methylation track showing methylation levels from BAM
pub fn draw_methylation_track(
    painter: &Painter,
    rect: Rect,
    track: &AlignmentTrack,
    viewport: &Viewport,
    y_offset: f32,
    height: f32,
) -> f32
{
    let track_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));

    // Background
    painter.rect_filled(track_rect, 0.0, Color32::from_gray(250));

    // Draw label
    painter.text(
        Pos2::new(rect.left() + 5.0, y_offset + 5.0),
        egui::Align2::LEFT_TOP,
        "Methylation",
        FontId::proportional(12.0),
        Color32::BLACK,
    );

    // Check if we have methylation data
    if track.methylation.is_empty()
    {
        painter.text(
            Pos2::new(rect.center().x, y_offset + height / 2.0),
            egui::Align2::CENTER_CENTER,
            "No methylation data",
            FontId::proportional(11.0),
            Color32::from_gray(150),
        );
        return track_rect.bottom();
    }

    // Get methylation data for viewport
    let bin_size = 100;
    let (region_start, region_end) = track.loaded_region.unwrap_or((0, track.reference_length));

    // Only show methylation if viewport overlaps with loaded region
    if viewport.end <= region_start || viewport.start >= region_end
    {
        painter.text(
            Pos2::new(rect.center().x, y_offset + height / 2.0),
            egui::Align2::CENTER_CENTER,
            "Methylation not loaded for this region",
            FontId::proportional(11.0),
            Color32::from_gray(150),
        );
        return track_rect.bottom();
    }

    // Calculate bin indices relative to the loaded region start
    let viewport_start_in_region = viewport.start.saturating_sub(region_start);
    let viewport_end_in_region = viewport.end.saturating_sub(region_start);

    let start_bin = viewport_start_in_region / bin_size;
    let end_bin = ((viewport_end_in_region + bin_size - 1) / bin_size).min(track.methylation.len());

    if start_bin >= track.methylation.len()
    {
        return track_rect.bottom();
    }

    let methylation_slice = &track.methylation[start_bin..end_bin];

    // Find max calls for scaling
    let max_calls = methylation_slice
        .iter()
        .map(|p| p.total_calls)
        .max()
        .unwrap_or(1)
        .max(1);

    // Count total methylation calls for label
    let total_meth_calls: u32 = methylation_slice.iter().map(|p| p.methylated_calls).sum();
    let total_calls: u32 = methylation_slice.iter().map(|p| p.total_calls).sum();

    // Update label with methylation info
    let avg_meth = if total_calls > 0
    {
        (total_meth_calls as f32 / total_calls as f32) * 100.0
    }
    else
    {
        0.0
    };

    painter.text(
        Pos2::new(rect.left() + 5.0, y_offset + 5.0),
        egui::Align2::LEFT_TOP,
        format!("Methylation ({:.1}% avg)", avg_meth),
        FontId::proportional(12.0),
        Color32::BLACK,
    );

    // Draw methylation bars
    let plot_height = height - 25.0;
    let y_base = track_rect.bottom() - 5.0;

    for point in methylation_slice
    {
        if point.position < viewport.start
            || point.position >= viewport.end
            || point.total_calls == 0
        {
            continue;
        }

        let x = viewport.position_to_screen(point.position, rect.width()) + rect.left();
        let bar_width = (rect.width() / viewport.width() as f32 * bin_size as f32).max(1.0);

        // Calculate methylation percentage
        let meth_pct = point.methylated_calls as f32 / point.total_calls as f32;

        // Height based on coverage (how many calls)
        let bar_height = (point.total_calls as f32 / max_calls as f32) * plot_height;

        // Color based on methylation percentage
        // Low methylation: blue, High methylation: red
        let r = (meth_pct * 255.0) as u8;
        let b = ((1.0 - meth_pct) * 255.0) as u8;
        let color = Color32::from_rgb(r, 50, b);

        painter.rect_filled(
            Rect::from_min_size(
                Pos2::new(x, y_base - bar_height),
                Vec2::new(bar_width, bar_height),
            ),
            0.0,
            color,
        );
    }

    // Draw legend
    let legend_y = y_offset + 18.0;
    let legend_x = rect.right() - 120.0;

    // Draw gradient legend
    for i in 0..50
    {
        let pct = i as f32 / 50.0;
        let r = (pct * 255.0) as u8;
        let b = ((1.0 - pct) * 255.0) as u8;
        painter.rect_filled(
            Rect::from_min_size(
                Pos2::new(legend_x + i as f32 * 2.0, legend_y),
                Vec2::new(2.0, 8.0),
            ),
            0.0,
            Color32::from_rgb(r, 50, b),
        );
    }

    painter.text(
        Pos2::new(legend_x - 5.0, legend_y + 4.0),
        egui::Align2::RIGHT_CENTER,
        "0%",
        FontId::proportional(9.0),
        Color32::DARK_GRAY,
    );

    painter.text(
        Pos2::new(legend_x + 105.0, legend_y + 4.0),
        egui::Align2::LEFT_CENTER,
        "100%",
        FontId::proportional(9.0),
        Color32::DARK_GRAY,
    );

    track_rect.bottom()
}

/// Draw TSV custom data track
pub fn draw_tsv_track(
    painter: &Painter,
    rect: Rect,
    track: &TsvChromosomeTrack,
    viewport: &Viewport,
    y_offset: f32,
    height: f32,
) -> f32
{
    let track_rect =
        Rect::from_min_size(Pos2::new(rect.left(), y_offset), Vec2::new(rect.width(), height));

    // Draw track background
    painter.rect_filled(track_rect, 0.0, Color32::from_gray(245));

    // Draw track border
    painter.rect_stroke(track_rect, 0.0, Stroke::new(1.0, Color32::from_gray(200)));

    // Draw track label
    painter.text(
        Pos2::new(rect.left() + 5.0, y_offset + 5.0),
        egui::Align2::LEFT_TOP,
        "Custom Track",
        FontId::proportional(12.0),
        Color32::from_gray(100),
    );

    // Get points in viewport range
    let points = track.get_points_in_region(viewport.start, viewport.end);

    if points.is_empty()
    {
        painter.text(
            Pos2::new(rect.center().x, y_offset + height / 2.0),
            egui::Align2::CENTER_CENTER,
            "No data in this region",
            FontId::proportional(12.0),
            Color32::from_gray(150),
        );
        return track_rect.bottom();
    }

    // Calculate value range for this view
    let min_signal = points
        .iter()
        .map(|p| p.signal)
        .fold(f64::INFINITY, f64::min);
    let max_signal = points
        .iter()
        .map(|p| p.signal)
        .fold(f64::NEG_INFINITY, f64::max);
    let signal_range = max_signal - min_signal;

    // Draw scale labels
    let label_y_top = y_offset + 20.0;
    let label_y_bottom = y_offset + height - 5.0;

    painter.text(
        Pos2::new(rect.left() + 5.0, label_y_top),
        egui::Align2::LEFT_TOP,
        format!("{:.2}", max_signal),
        FontId::proportional(10.0),
        Color32::from_gray(120),
    );

    painter.text(
        Pos2::new(rect.left() + 5.0, label_y_bottom),
        egui::Align2::LEFT_BOTTOM,
        format!("{:.2}", min_signal),
        FontId::proportional(10.0),
        Color32::from_gray(120),
    );

    // Draw zero line if range includes zero
    if min_signal <= 0.0 && max_signal >= 0.0 && signal_range > 0.0
    {
        let zero_y =
            y_offset + 25.0 + (height - 30.0) * (1.0 - ((0.0 - min_signal) / signal_range) as f32);
        painter.line_segment(
            [
                Pos2::new(rect.left() + 50.0, zero_y),
                Pos2::new(rect.right(), zero_y),
            ],
            Stroke::new(1.0, Color32::from_rgb(150, 150, 150)),
        );
    }

    // Draw data as line graph
    let data_area_top = y_offset + 25.0;
    let data_area_height = height - 30.0;

    let mut prev_pos: Option<Pos2> = None;

    for point in &points
    {
        let x = rect.left()
            + ((point.position - viewport.start) as f32 / (viewport.end - viewport.start) as f32)
                * rect.width();

        let normalized_signal = if signal_range > 0.0
        {
            ((point.signal - min_signal) / signal_range) as f32
        }
        else
        {
            0.5
        };

        let y = data_area_top + data_area_height * (1.0 - normalized_signal);
        let pos = Pos2::new(x, y);

        // Draw point
        painter.circle_filled(pos, 2.0, Color32::from_rgb(50, 100, 200));

        // Draw line to previous point
        if let Some(prev) = prev_pos
        {
            painter.line_segment([prev, pos], Stroke::new(1.5, Color32::from_rgb(50, 100, 200)));
        }

        prev_pos = Some(pos);
    }

    track_rect.bottom()
}
