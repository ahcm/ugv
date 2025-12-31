mod fasta;
mod gff;
mod interval_tree;
mod renderer;
mod viewport;

use anyhow::Result;
use eframe::egui;
use std::path::PathBuf;

fn main() -> Result<()>
{
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_title("Ultra-Fast Genome Viewer"),
        ..Default::default()
    };

    eframe::run_native("UGV", options, Box::new(|_cc| Ok(Box::new(GenomeViewer::new()))))
        .map_err(|e| anyhow::anyhow!("Failed to run app: {}", e))
}

struct GenomeViewer
{
    genome: Option<fasta::Genome>,
    features: Vec<gff::Feature>,
    interval_tree: Option<interval_tree::IntervalTree>,
    viewport: viewport::Viewport,
    fasta_path: String,
    gff_path: String,
    selected_chromosome: Option<String>,
    status_message: String,
}

impl GenomeViewer
{
    fn new() -> Self
    {
        Self {
            genome: None,
            features: Vec::new(),
            interval_tree: None,
            viewport: viewport::Viewport::new(0, 1000000),
            fasta_path: String::new(),
            gff_path: String::new(),
            selected_chromosome: None,
            status_message: "Load a FASTA file to begin".to_string(),
        }
    }

    fn load_fasta(&mut self)
    {
        if self.fasta_path.is_empty()
        {
            self.status_message = "Error: No FASTA path specified".to_string();
            return;
        }

        match fasta::Genome::from_file(&self.fasta_path)
        {
            Ok(genome) =>
            {
                self.status_message = format!("Loaded {} chromosomes", genome.chromosomes.len());
                if let Some(first_chr) = genome.chromosomes.keys().next()
                {
                    self.selected_chromosome = Some(first_chr.clone());
                    if let Some(chr) = genome.chromosomes.get(first_chr)
                    {
                        self.viewport = viewport::Viewport::new(0, chr.length);
                    }
                }
                self.genome = Some(genome);
            }
            Err(e) =>
            {
                self.status_message = format!("Error loading FASTA: {}", e);
            }
        }
    }

    fn load_gff(&mut self)
    {
        if self.gff_path.is_empty()
        {
            self.status_message = "Error: No GFF path specified".to_string();
            return;
        }

        match gff::parse_gff(&self.gff_path)
        {
            Ok(features) =>
            {
                let count = features.len();
                self.interval_tree = Some(interval_tree::IntervalTree::from_features(&features));
                self.features = features;
                self.status_message = format!("Loaded {} features", count);
            }
            Err(e) =>
            {
                self.status_message = format!("Error loading GFF: {}", e);
            }
        }
    }
}

impl eframe::App for GenomeViewer
{
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame)
    {
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("Genome Viewer");
                ui.separator();
                ui.label("FASTA:");
                ui.text_edit_singleline(&mut self.fasta_path);
                if ui.button("Load FASTA").clicked()
                {
                    self.load_fasta();
                }
                ui.separator();
                ui.label("GFF:");
                ui.text_edit_singleline(&mut self.gff_path);
                if ui.button("Load GFF").clicked()
                {
                    self.load_gff();
                }
            });
        });

        egui::TopBottomPanel::bottom("bottom_panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(&self.status_message);
                ui.separator();
                if let Some(ref chr) = self.selected_chromosome
                {
                    ui.label(format!(
                        "Chromosome: {} | Position: {}-{}",
                        chr, self.viewport.start, self.viewport.end
                    ));
                }
            });
        });

        egui::SidePanel::left("side_panel")
            .min_width(200.0)
            .show(ctx, |ui| {
                ui.heading("Chromosomes");
                ui.separator();

                if let Some(ref genome) = self.genome
                {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        for (chr_name, chr_data) in &genome.chromosomes
                        {
                            let is_selected = self.selected_chromosome.as_ref() == Some(chr_name);
                            if ui
                                .selectable_label(
                                    is_selected,
                                    format!("{} ({} bp)", chr_name, chr_data.length),
                                )
                                .clicked()
                            {
                                self.selected_chromosome = Some(chr_name.clone());
                                self.viewport = viewport::Viewport::new(0, chr_data.length);
                            }
                        }
                    });
                }
                else
                {
                    ui.label("No genome loaded");
                }
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            let available_size = ui.available_size();
            let (response, painter) =
                ui.allocate_painter(available_size, egui::Sense::click_and_drag());

            // Handle pan and zoom
            if response.dragged()
            {
                let delta = response.drag_delta();
                let pixels_per_base =
                    available_size.x / (self.viewport.end - self.viewport.start) as f32;
                let bases_moved = (delta.x / pixels_per_base) as i64;
                self.viewport.pan(-bases_moved);
            }

            if let Some(hover_pos) = response.hover_pos()
            {
                let scroll_delta = ui.input(|i| i.smooth_scroll_delta.y);
                if scroll_delta != 0.0
                {
                    let mouse_x = hover_pos.x - response.rect.left();
                    let relative_pos = mouse_x / available_size.x;
                    let zoom_factor = if scroll_delta > 0.0 { 0.9 } else { 1.1 };
                    self.viewport.zoom(zoom_factor, relative_pos);
                }
            }

            // Render genome
            if let (Some(ref genome), Some(ref chr_name)) =
                (&self.genome, &self.selected_chromosome)
            {
                if let Some(chr) = genome.chromosomes.get(chr_name)
                {
                    renderer::render_genome(
                        &painter,
                        response.rect,
                        chr,
                        &self.viewport,
                        &self.features,
                        self.interval_tree.as_ref(),
                        chr_name,
                    );
                }
            }
            else
            {
                painter.text(
                    response.rect.center(),
                    egui::Align2::CENTER_CENTER,
                    "Load a genome file to begin",
                    egui::FontId::proportional(20.0),
                    egui::Color32::GRAY,
                );
            }
        });

        ctx.request_repaint();
    }
}
