use crate::fasta::Chromosome;
use crate::gff::Feature;
use crate::interval_tree::IntervalTree;
use crate::viewport::Viewport;
use egui::{Color32, Painter, Pos2, Rect, Stroke, Vec2};

const RULER_HEIGHT: f32 = 30.0;
const SEQUENCE_HEIGHT: f32 = 40.0;
const GC_CONTENT_HEIGHT: f32 = 60.0;
const FEATURE_TRACK_HEIGHT: f32 = 30.0;
const TRACK_SPACING: f32 = 10.0;

pub fn render_genome(
    painter: &Painter,
    rect: Rect,
    chromosome: &Chromosome,
    viewport: &Viewport,
    features: &[Feature],
    interval_tree: Option<&IntervalTree>,
    chr_name: &str,
)
{
    let mut y_offset = rect.top();

    // Draw ruler
    y_offset = draw_ruler(painter, rect, viewport, y_offset);
    y_offset += TRACK_SPACING;

    // Draw GC content plot
    y_offset = draw_gc_content(painter, rect, chromosome, viewport, y_offset);
    y_offset += TRACK_SPACING;

    // Draw sequence (if zoomed in enough)
    if viewport.width() < 1000
    {
        y_offset = draw_sequence(painter, rect, chromosome, viewport, y_offset);
        y_offset += TRACK_SPACING;
    }

    // Draw features
    if let Some(tree) = interval_tree
    {
        draw_features(painter, rect, features, viewport, tree, chr_name, y_offset);
    }
}

fn draw_ruler(painter: &Painter, rect: Rect, viewport: &Viewport, y_offset: f32) -> f32
{
    let ruler_rect = Rect::from_min_size(
        Pos2::new(rect.left(), y_offset),
        Vec2::new(rect.width(), RULER_HEIGHT),
    );

    // Background
    painter.rect_filled(ruler_rect, 0.0, Color32::from_gray(240));

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
        let label = format_position(tick);
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

fn draw_gc_content(
    painter: &Painter,
    rect: Rect,
    chromosome: &Chromosome,
    viewport: &Viewport,
    y_offset: f32,
) -> f32
{
    let gc_rect = Rect::from_min_size(
        Pos2::new(rect.left(), y_offset),
        Vec2::new(rect.width(), GC_CONTENT_HEIGHT),
    );

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
        let y = gc_rect.bottom() - (gc * GC_CONTENT_HEIGHT);
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

fn draw_sequence(
    painter: &Painter,
    rect: Rect,
    chromosome: &Chromosome,
    viewport: &Viewport,
    y_offset: f32,
) -> f32
{
    let seq_rect = Rect::from_min_size(
        Pos2::new(rect.left(), y_offset),
        Vec2::new(rect.width(), SEQUENCE_HEIGHT),
    );

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

fn draw_features(
    painter: &Painter,
    rect: Rect,
    features: &[Feature],
    viewport: &Viewport,
    interval_tree: &IntervalTree,
    chr_name: &str,
    y_offset: f32,
)
{
    // Query features in viewport
    let feature_indices = interval_tree.query(chr_name, viewport.start, viewport.end);

    if feature_indices.is_empty()
    {
        return;
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

    // Draw each track
    for (track_idx, track) in tracks.iter().enumerate()
    {
        let track_y = y_offset + (track_idx as f32 * (FEATURE_TRACK_HEIGHT + TRACK_SPACING));

        for &idx in track
        {
            let feature = &features[idx];
            draw_feature(painter, rect, feature, viewport, track_y);
        }
    }
}

fn draw_feature(
    painter: &Painter,
    rect: Rect,
    feature: &Feature,
    viewport: &Viewport,
    y_offset: f32,
)
{
    let start_x = viewport.position_to_screen(feature.start, rect.width()) + rect.left();
    let end_x = viewport.position_to_screen(feature.end, rect.width()) + rect.left();

    let feature_rect = Rect::from_min_size(
        Pos2::new(start_x, y_offset),
        Vec2::new((end_x - start_x).max(2.0), FEATURE_TRACK_HEIGHT),
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

fn format_position(pos: usize) -> String
{
    if pos >= 1_000_000
    {
        format!("{:.1}M", pos as f32 / 1_000_000.0)
    }
    else if pos >= 1_000
    {
        format!("{:.1}K", pos as f32 / 1_000.0)
    }
    else
    {
        pos.to_string()
    }
}
