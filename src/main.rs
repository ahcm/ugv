mod bam;
mod fasta;
mod file_loader;
mod gff;
mod interval_tree;
mod renderer;
mod translation;
mod viewport;

use anyhow::Result;
use eframe::egui;

#[cfg(not(target_arch = "wasm32"))]
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

#[cfg(target_arch = "wasm32")]
fn main()
{
    use eframe::wasm_bindgen::JsCast as _;

    // Redirect tracing to console.log
    eframe::WebLogger::init(log::LevelFilter::Debug).ok();

    let web_options = eframe::WebOptions::default();

    wasm_bindgen_futures::spawn_local(async {
        let document = web_sys::window()
            .expect("No window")
            .document()
            .expect("No document");

        let canvas = document
            .get_element_by_id("the_canvas_id")
            .expect("Failed to find the_canvas_id")
            .dyn_into::<web_sys::HtmlCanvasElement>()
            .expect("the_canvas_id was not a HtmlCanvasElement");

        eframe::WebRunner::new()
            .start(canvas, web_options, Box::new(|_cc| Ok(Box::new(GenomeViewer::new()))))
            .await
            .expect("failed to start eframe");
    });
}

#[derive(PartialEq, Clone, Copy)]
enum ChromosomeSort
{
    Natural,
    Alphabetical,
    Size,
}

type GenomePromise = egui::mutex::Mutex<Option<poll_promise::Promise<Result<fasta::Genome>>>>;
type FeaturesPromise = egui::mutex::Mutex<Option<poll_promise::Promise<Result<Vec<gff::Feature>>>>>;
type BamPromise = egui::mutex::Mutex<Option<poll_promise::Promise<Result<bam::AlignmentData>>>>;

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
    chromosome_search: String,
    chromosome_sort: ChromosomeSort,
    position_search: String,
    feature_search: String,
    search_results: Vec<(String, String, usize, usize)>, // (chr, name, start, end)
    show_search_results: bool,
    show_amino_acids: bool,
    // BAM support
    alignments: Option<bam::AlignmentData>,
    bam_path: String,
    show_coverage: bool,
    show_alignments: bool,
    show_variants: bool,
    #[cfg(target_arch = "wasm32")]
    genome_promise: GenomePromise,
    #[cfg(target_arch = "wasm32")]
    features_promise: FeaturesPromise,
    #[cfg(target_arch = "wasm32")]
    fasta_file_name: String,
    #[cfg(target_arch = "wasm32")]
    gff_file_name: String,
    #[cfg(target_arch = "wasm32")]
    bam_promise: BamPromise,
    #[cfg(target_arch = "wasm32")]
    bam_file_name: String,
    #[cfg(target_arch = "wasm32")]
    file_picker_open: Option<FilePickerType>,
    loading_progress: Option<LoadingProgress>,
}

#[cfg(target_arch = "wasm32")]
#[derive(Clone, Copy, PartialEq)]
enum FilePickerType
{
    Fasta,
    Gff,
    Bam,
}

#[derive(Clone)]
struct LoadingProgress
{
    description: String,
    bytes_loaded: usize,
    total_bytes: Option<usize>,
}

fn format_bytes(bytes: usize) -> String
{
    const KB: usize = 1024;
    const MB: usize = KB * 1024;
    const GB: usize = MB * 1024;

    if bytes >= GB
    {
        format!("{:.2} GB", bytes as f64 / GB as f64)
    }
    else if bytes >= MB
    {
        format!("{:.2} MB", bytes as f64 / MB as f64)
    }
    else if bytes >= KB
    {
        format!("{:.2} KB", bytes as f64 / KB as f64)
    }
    else
    {
        format!("{} B", bytes)
    }
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
            chromosome_search: String::new(),
            chromosome_sort: ChromosomeSort::Natural,
            position_search: String::new(),
            feature_search: String::new(),
            search_results: Vec::new(),
            show_search_results: false,
            show_amino_acids: false,
            // BAM support
            alignments: None,
            bam_path: String::new(),
            show_coverage: true,
            show_alignments: true,
            show_variants: false,
            #[cfg(target_arch = "wasm32")]
            genome_promise: egui::mutex::Mutex::new(None),
            #[cfg(target_arch = "wasm32")]
            features_promise: egui::mutex::Mutex::new(None),
            #[cfg(target_arch = "wasm32")]
            fasta_file_name: String::new(),
            #[cfg(target_arch = "wasm32")]
            gff_file_name: String::new(),
            #[cfg(target_arch = "wasm32")]
            bam_promise: egui::mutex::Mutex::new(None),
            #[cfg(target_arch = "wasm32")]
            bam_file_name: String::new(),
            #[cfg(target_arch = "wasm32")]
            file_picker_open: None,
            loading_progress: None,
        }
    }

    fn extract_chromosome_number(name: &str) -> Option<u32>
    {
        // Try to extract numeric part from chromosome names like "chr1", "1", "chromosome1"
        let s = name.to_lowercase();
        let s = s.strip_prefix("chr").unwrap_or(&s);
        let s = s.strip_prefix("chromosome").unwrap_or(s);

        // Extract leading digits
        let digits: String = s.chars().take_while(|c| c.is_ascii_digit()).collect();
        digits.parse().ok()
    }

    fn chromosome_natural_sort_key(name: &str) -> (u8, u32, String)
    {
        // Sort order: autosomes (1-22), sex chromosomes (X, Y), mitochondrial (M, MT), others
        let lower = name.to_lowercase();

        if let Some(num) = Self::extract_chromosome_number(name)
        {
            // Numeric chromosomes come first
            (0, num, String::new())
        }
        else if lower.contains("x") || lower == "chrx"
        {
            (1, 0, String::from("x"))
        }
        else if lower.contains("y") || lower == "chry"
        {
            (1, 1, String::from("y"))
        }
        else if lower.contains("m") || lower == "chrm" || lower == "mt" || lower == "chrmt"
        {
            (1, 2, String::from("m"))
        }
        else
        {
            // Everything else sorted alphabetically at the end
            (2, 0, name.to_lowercase())
        }
    }

    fn parse_position_search(&self, query: &str) -> Option<(String, usize)>
    {
        let query = query.trim();

        // Format: "chr1:100000" or "chr1:100,000" or just "100000"
        if let Some((chr, pos)) = query.split_once(':')
        {
            // Has chromosome specified
            let pos_str = pos.replace(&[',', '_'][..], "");
            if let Ok(position) = pos_str.parse::<usize>()
            {
                return Some((chr.to_string(), position));
            }
        }
        else
        {
            // Just a number, use current chromosome
            let pos_str = query.replace(&[',', '_'][..], "");
            if let Ok(position) = pos_str.parse::<usize>()
            {
                if let Some(ref chr) = self.selected_chromosome
                {
                    return Some((chr.clone(), position));
                }
            }
        }

        None
    }

    fn search_position(&mut self)
    {
        if let Some((chr, position)) = self.parse_position_search(&self.position_search)
        {
            // Check if chromosome exists
            if let Some(ref genome) = self.genome
            {
                if let Some(chr_data) = genome.chromosomes.get(&chr)
                {
                    self.selected_chromosome = Some(chr.clone());

                    // Center the view on the position
                    let window_size = 10000.min(chr_data.length / 2);
                    let start = position.saturating_sub(window_size / 2);
                    let end = (position + window_size / 2).min(chr_data.length);

                    self.viewport = viewport::Viewport::new(0, chr_data.length);
                    self.viewport.set_region(start, end);

                    self.status_message = format!("Jumped to {}:{}", chr, position);
                }
                else
                {
                    self.status_message = format!("Chromosome '{}' not found", chr);
                }
            }
        }
        else
        {
            self.status_message =
                "Invalid position format. Use 'chr1:100000' or '100000'".to_string();
        }
    }

    fn search_features(&mut self)
    {
        let query = self.feature_search.trim().to_lowercase();
        if query.is_empty()
        {
            return;
        }

        self.search_results.clear();

        for feature in &self.features
        {
            let name = feature.name().to_lowercase();
            let id = feature.get_attribute("ID").unwrap_or("").to_lowercase();
            let gene = feature.get_attribute("gene").unwrap_or("").to_lowercase();

            if name.contains(&query) || id.contains(&query) || gene.contains(&query)
            {
                self.search_results.push((
                    feature.seqid.clone(),
                    feature.name(),
                    feature.start,
                    feature.end,
                ));
            }
        }

        if self.search_results.is_empty()
        {
            self.status_message = format!("No features found matching '{}'", self.feature_search);
        }
        else
        {
            self.status_message = format!(
                "Found {} feature(s) matching '{}'",
                self.search_results.len(),
                self.feature_search
            );
            self.show_search_results = true;
        }
    }

    fn jump_to_feature(&mut self, chr: &str, start: usize, end: usize)
    {
        if let Some(ref genome) = self.genome
        {
            if let Some(chr_data) = genome.chromosomes.get(chr)
            {
                self.selected_chromosome = Some(chr.to_string());

                // Zoom to feature with some padding
                self.viewport = viewport::Viewport::new(0, chr_data.length);
                self.viewport.zoom_to_feature(start, end, 0.2);

                self.show_search_results = false;
                self.status_message = format!("Jumped to {}:{}-{}", chr, start, end);
            }
        }
    }

    fn get_sorted_chromosomes(&self) -> Vec<(String, usize)>
    {
        if let Some(ref genome) = self.genome
        {
            let mut chromosomes: Vec<_> = genome
                .chromosomes
                .iter()
                .map(|(name, chr)| (name.clone(), chr.length))
                .collect();

            // Apply search filter
            if !self.chromosome_search.is_empty()
            {
                let search_lower = self.chromosome_search.to_lowercase();
                chromosomes.retain(|(name, _)| name.to_lowercase().contains(&search_lower));
            }

            // Apply sort
            match self.chromosome_sort
            {
                ChromosomeSort::Natural =>
                {
                    chromosomes.sort_by(|a, b| {
                        let key_a = Self::chromosome_natural_sort_key(&a.0);
                        let key_b = Self::chromosome_natural_sort_key(&b.0);
                        key_a.cmp(&key_b)
                    });
                }
                ChromosomeSort::Alphabetical =>
                {
                    chromosomes.sort_by(|a, b| a.0.to_lowercase().cmp(&b.0.to_lowercase()));
                }
                ChromosomeSort::Size =>
                {
                    chromosomes.sort_by(|a, b| b.1.cmp(&a.1));
                }
            }

            chromosomes
        }
        else
        {
            Vec::new()
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
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

    #[cfg(target_arch = "wasm32")]
    fn load_fasta(&mut self)
    {
        if self.fasta_path.is_empty()
        {
            self.status_message = "Error: No FASTA path/URL specified".to_string();
            return;
        }

        let path = self.fasta_path.clone();
        self.status_message = format!("Loading FASTA from {}...", path);
        self.loading_progress = Some(LoadingProgress {
            description: "Downloading FASTA file...".to_string(),
            bytes_loaded: 0,
            total_bytes: None,
        });

        let promise = poll_promise::Promise::spawn_local(async move {
            let data = file_loader::load_file_async(&path).await?;
            let is_gzipped = path.ends_with(".gz");
            fasta::Genome::from_bytes(data, is_gzipped)
        });

        *self.genome_promise.lock() = Some(promise);
    }

    #[cfg(target_arch = "wasm32")]
    fn load_fasta_from_file(&mut self, file: web_sys::File)
    {
        let file_name = file.name();
        let file_size = file.size() as usize;
        self.fasta_file_name = file_name.clone();
        self.status_message = format!("Loading FASTA from {}...", file_name);
        self.loading_progress = Some(LoadingProgress {
            description: "Reading FASTA file...".to_string(),
            bytes_loaded: 0,
            total_bytes: Some(file_size),
        });

        let is_gzipped = file_name.ends_with(".gz");

        let promise = poll_promise::Promise::spawn_local(async move {
            use gloo_file::futures::read_as_bytes;

            let data = read_as_bytes(&file.into())
                .await
                .map_err(|e| anyhow::anyhow!("Failed to read file: {:?}", e))?;

            fasta::Genome::from_bytes(data, is_gzipped)
        });

        *self.genome_promise.lock() = Some(promise);
    }

    #[cfg(target_arch = "wasm32")]
    fn check_genome_promise(&mut self)
    {
        let mut promise_lock = self.genome_promise.lock();
        if let Some(promise) = promise_lock.as_mut()
        {
            if let Some(result) = promise.ready()
            {
                self.loading_progress = None;
                match result
                {
                    Ok(genome) =>
                    {
                        self.status_message =
                            format!("Loaded {} chromosomes", genome.chromosomes.len());
                        if let Some(first_chr) = genome.chromosomes.keys().next()
                        {
                            self.selected_chromosome = Some(first_chr.clone());
                            if let Some(chr) = genome.chromosomes.get(first_chr)
                            {
                                self.viewport = viewport::Viewport::new(0, chr.length);
                            }
                        }
                        self.genome = Some(genome.clone());
                    }
                    Err(e) =>
                    {
                        self.status_message = format!("Error loading FASTA: {}", e);
                    }
                }
                *promise_lock = None;
            }
            else
            {
                // Still loading - update progress description
                if let Some(ref mut progress) = self.loading_progress
                {
                    // Update description based on how long we've been loading
                    progress.description = "Parsing FASTA...".to_string();
                }
            }
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
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

    #[cfg(target_arch = "wasm32")]
    fn load_gff(&mut self)
    {
        if self.gff_path.is_empty()
        {
            self.status_message = "Error: No GFF path/URL specified".to_string();
            return;
        }

        let path = self.gff_path.clone();
        self.status_message = format!("Loading GFF from {}...", path);
        self.loading_progress = Some(LoadingProgress {
            description: "Downloading GFF file...".to_string(),
            bytes_loaded: 0,
            total_bytes: None,
        });

        let promise = poll_promise::Promise::spawn_local(async move {
            let data = file_loader::load_file_async(&path).await?;
            let is_gzipped = path.ends_with(".gz");
            gff::parse_gff_bytes(data, is_gzipped)
        });

        *self.features_promise.lock() = Some(promise);
    }

    #[cfg(target_arch = "wasm32")]
    fn load_gff_from_file(&mut self, file: web_sys::File)
    {
        let file_name = file.name();
        let file_size = file.size() as usize;
        self.gff_file_name = file_name.clone();
        self.status_message = format!("Loading GFF from {}...", file_name);
        self.loading_progress = Some(LoadingProgress {
            description: "Reading GFF file...".to_string(),
            bytes_loaded: 0,
            total_bytes: Some(file_size),
        });

        let is_gzipped = file_name.ends_with(".gz");

        let promise = poll_promise::Promise::spawn_local(async move {
            use gloo_file::futures::read_as_bytes;

            let data = read_as_bytes(&file.into())
                .await
                .map_err(|e| anyhow::anyhow!("Failed to read file: {:?}", e))?;

            gff::parse_gff_bytes(data, is_gzipped)
        });

        *self.features_promise.lock() = Some(promise);
    }

    // BAM loading functions
    #[cfg(not(target_arch = "wasm32"))]
    fn load_bam(&mut self)
    {
        if self.bam_path.is_empty()
        {
            self.status_message = "Error: No BAM path specified".to_string();
            return;
        }

        match bam::AlignmentData::from_file(&self.bam_path)
        {
            Ok(alignments) =>
            {
                let total_reads: usize = alignments.tracks.values().map(|t| t.records.len()).sum();
                self.status_message = format!(
                    "Loaded {} reads from {} references",
                    total_reads,
                    alignments.tracks.len()
                );
                self.alignments = Some(alignments);
            }
            Err(e) =>
            {
                self.status_message = format!("Error loading BAM: {}", e);
            }
        }
    }

    #[cfg(target_arch = "wasm32")]
    fn load_bam(&mut self)
    {
        if self.bam_path.is_empty()
        {
            self.status_message = "Error: No BAM path/URL specified".to_string();
            return;
        }

        let path = self.bam_path.clone();
        self.status_message = format!("Loading BAM from {}...", path);
        self.loading_progress = Some(LoadingProgress {
            description: "Downloading BAM file...".to_string(),
            bytes_loaded: 0,
            total_bytes: None,
        });

        let promise = poll_promise::Promise::spawn_local(async move {
            let data = file_loader::load_file_async(&path).await?;
            bam::AlignmentData::from_bytes(data)
        });

        *self.bam_promise.lock() = Some(promise);
    }

    #[cfg(target_arch = "wasm32")]
    fn load_bam_from_file(&mut self, file: web_sys::File)
    {
        let file_name = file.name();
        let file_size = file.size() as usize;
        self.bam_file_name = file_name.clone();
        self.status_message = format!("Loading BAM from {}...", file_name);
        self.loading_progress = Some(LoadingProgress {
            description: "Reading BAM file...".to_string(),
            bytes_loaded: 0,
            total_bytes: Some(file_size),
        });

        let promise = poll_promise::Promise::spawn_local(async move {
            use gloo_file::futures::read_as_bytes;

            let data = read_as_bytes(&file.into())
                .await
                .map_err(|e| anyhow::anyhow!("Failed to read file: {:?}", e))?;

            bam::AlignmentData::from_bytes(data)
        });

        *self.bam_promise.lock() = Some(promise);
    }

    #[cfg(target_arch = "wasm32")]
    fn open_file_picker(&mut self, picker_type: FilePickerType)
    {
        use wasm_bindgen::JsCast;
        use wasm_bindgen::closure::Closure;

        let window = web_sys::window().expect("no global window");
        let document = window.document().expect("no document");

        // Create file input element
        let input = document
            .create_element("input")
            .expect("failed to create input")
            .dyn_into::<web_sys::HtmlInputElement>()
            .expect("not an input element");

        input.set_type("file");

        // Set accept attribute based on file type
        match picker_type
        {
            FilePickerType::Fasta =>
            {
                input.set_accept(".fasta,.fa,.fna,.ffn,.faa,.frn,.fasta.gz,.fa.gz,.fna.gz,.ffn.gz,.faa.gz,.frn.gz,.gz");
            }
            FilePickerType::Gff =>
            {
                input.set_accept(".gff,.gff3,.gtf,.gff.gz,.gff3.gz,.gtf.gz,.gz");
            }
            FilePickerType::Bam =>
            {
                input.set_accept(".bam");
            }
        }

        // Store the picker type so we know what to do with the selected file
        self.file_picker_open = Some(picker_type);

        // Create a change event handler
        let input_clone = input.clone();
        let closure = Closure::wrap(Box::new(move |_event: web_sys::Event| {
            if let Some(files) = input_clone.files()
            {
                if files.length() > 0
                {
                    if let Some(file) = files.get(0)
                    {
                        // File selected - we'll handle this in the update loop
                        // by checking input.files()
                        log::info!("File selected: {}", file.name());
                    }
                }
            }
        }) as Box<dyn FnMut(_)>);

        input.set_onchange(Some(closure.as_ref().unchecked_ref()));
        closure.forget(); // Keep the closure alive

        // Trigger the file picker
        input.click();

        // Store the input element in the DOM so we can check it later
        // We'll attach it to the body but make it hidden
        input.style().set_property("display", "none").ok();
        input.set_id(&format!(
            "file_picker_{}",
            match picker_type
            {
                FilePickerType::Fasta => "fasta",
                FilePickerType::Gff => "gff",
                FilePickerType::Bam => "bam",
            }
        ));

        document
            .body()
            .expect("no body")
            .append_child(&input)
            .expect("failed to append input");
    }

    #[cfg(target_arch = "wasm32")]
    fn check_file_picker(&mut self)
    {
        use wasm_bindgen::JsCast;

        if let Some(picker_type) = self.file_picker_open
        {
            let window = web_sys::window().expect("no global window");
            let document = window.document().expect("no document");

            let element_id = format!(
                "file_picker_{}",
                match picker_type
                {
                    FilePickerType::Fasta => "fasta",
                    FilePickerType::Gff => "gff",
                    FilePickerType::Bam => "bam",
                }
            );

            if let Some(element) = document.get_element_by_id(&element_id)
            {
                let element_clone = element.clone();
                if let Ok(input) = element_clone.dyn_into::<web_sys::HtmlInputElement>()
                {
                    if let Some(files) = input.files()
                    {
                        if files.length() > 0
                        {
                            if let Some(file) = files.get(0)
                            {
                                // File was selected, load it
                                match picker_type
                                {
                                    FilePickerType::Fasta =>
                                    {
                                        self.load_fasta_from_file(file);
                                    }
                                    FilePickerType::Gff =>
                                    {
                                        self.load_gff_from_file(file);
                                    }
                                    FilePickerType::Bam =>
                                    {
                                        self.load_bam_from_file(file);
                                    }
                                }

                                // Remove the input element from DOM
                                if let Some(parent) = element.parent_element()
                                {
                                    parent.remove_child(&element).ok();
                                }

                                self.file_picker_open = None;
                            }
                        }
                    }
                }
            }
        }
    }

    #[cfg(target_arch = "wasm32")]
    fn check_features_promise(&mut self)
    {
        let mut promise_lock = self.features_promise.lock();
        if let Some(promise) = promise_lock.as_mut()
        {
            if let Some(result) = promise.ready()
            {
                self.loading_progress = None;
                match result
                {
                    Ok(features) =>
                    {
                        let count = features.len();
                        self.interval_tree =
                            Some(interval_tree::IntervalTree::from_features(&features));
                        self.features = features.clone();
                        self.status_message = format!("Loaded {} features", count);
                    }
                    Err(e) =>
                    {
                        self.status_message = format!("Error loading GFF: {}", e);
                    }
                }
                *promise_lock = None;
            }
            else
            {
                // Still loading - update progress description
                if let Some(ref mut progress) = self.loading_progress
                {
                    progress.description = "Parsing GFF...".to_string();
                }
            }
        }
    }

    #[cfg(target_arch = "wasm32")]
    fn check_bam_promise(&mut self)
    {
        let mut promise_lock = self.bam_promise.lock();
        if let Some(promise) = promise_lock.as_mut()
        {
            if let Some(result) = promise.ready()
            {
                self.loading_progress = None;
                match result
                {
                    Ok(alignments) =>
                    {
                        let total_reads: usize =
                            alignments.tracks.values().map(|t| t.records.len()).sum();
                        self.status_message = format!(
                            "Loaded {} reads from {} references",
                            total_reads,
                            alignments.tracks.len()
                        );
                        self.alignments = Some(alignments.clone());
                    }
                    Err(e) =>
                    {
                        self.status_message = format!("Error loading BAM: {}", e);
                    }
                }
                *promise_lock = None;
            }
            else
            {
                // Still loading - update progress description
                if let Some(ref mut progress) = self.loading_progress
                {
                    progress.description = "Parsing BAM...".to_string();
                }
            }
        }
    }
}

impl eframe::App for GenomeViewer
{
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame)
    {
        // Check for async file loading completion (WASM only)
        #[cfg(target_arch = "wasm32")]
        {
            self.check_file_picker();
            self.check_genome_promise();
            self.check_features_promise();
            self.check_bam_promise();

            // Handle drag and drop files
            ctx.input(|i| {
                if !i.raw.dropped_files.is_empty()
                {
                    for file in &i.raw.dropped_files
                    {
                        if let Some(bytes) = &file.bytes
                        {
                            // File already loaded as bytes
                            let name = file.name.clone();
                            let file_size = bytes.len();
                            let is_gzipped = name.ends_with(".gz");

                            // Detect file type by extension
                            let lower_name = name.to_lowercase();
                            if lower_name.contains(".fasta")
                                || lower_name.contains(".fa")
                                || lower_name.contains(".fna")
                                || lower_name.contains(".ffn")
                                || lower_name.contains(".faa")
                                || lower_name.contains(".frn")
                            {
                                // Load as FASTA
                                self.fasta_file_name = name.clone();
                                self.status_message = format!("Loading FASTA from {}...", name);
                                self.loading_progress = Some(LoadingProgress {
                                    description: "Parsing FASTA...".to_string(),
                                    bytes_loaded: file_size,
                                    total_bytes: Some(file_size),
                                });

                                let data = bytes.clone().to_vec();
                                let promise = poll_promise::Promise::spawn_local(async move {
                                    fasta::Genome::from_bytes(data, is_gzipped)
                                });
                                *self.genome_promise.lock() = Some(promise);
                            }
                            else if lower_name.contains(".gff") || lower_name.contains(".gtf")
                            {
                                // Load as GFF/GTF
                                self.gff_file_name = name.clone();
                                self.status_message = format!("Loading GFF from {}...", name);
                                self.loading_progress = Some(LoadingProgress {
                                    description: "Parsing GFF...".to_string(),
                                    bytes_loaded: file_size,
                                    total_bytes: Some(file_size),
                                });

                                let data = bytes.clone().to_vec();
                                let promise = poll_promise::Promise::spawn_local(async move {
                                    gff::parse_gff_bytes(data, is_gzipped)
                                });
                                *self.features_promise.lock() = Some(promise);
                            }
                            else if lower_name.ends_with(".bam")
                            {
                                // Load as BAM
                                self.bam_file_name = name.clone();
                                self.status_message = format!("Loading BAM from {}...", name);
                                self.loading_progress = Some(LoadingProgress {
                                    description: "Parsing BAM...".to_string(),
                                    bytes_loaded: file_size,
                                    total_bytes: Some(file_size),
                                });

                                let data = bytes.clone().to_vec();
                                let promise = poll_promise::Promise::spawn_local(async move {
                                    bam::AlignmentData::from_bytes(data)
                                });
                                *self.bam_promise.lock() = Some(promise);
                            }
                        }
                    }
                }
            });
        }

        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("Genome Viewer");
                ui.separator();

                #[cfg(not(target_arch = "wasm32"))]
                {
                    if ui.button("Open FASTA...").clicked()
                    {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("FASTA", &["fasta", "fa", "fna", "ffn", "faa", "frn"])
                            .add_filter(
                                "FASTA (gzipped)",
                                &[
                                    "fasta.gz", "fa.gz", "fna.gz", "ffn.gz", "faa.gz", "frn.gz",
                                    "gz",
                                ],
                            )
                            .add_filter("All files", &["*"])
                            .pick_file()
                        {
                            self.fasta_path = path.display().to_string();
                            self.load_fasta();
                        }
                    }

                    if !self.fasta_path.is_empty()
                    {
                        ui.label(format!(
                            "📄 {}",
                            std::path::Path::new(&self.fasta_path)
                                .file_name()
                                .and_then(|n| n.to_str())
                                .unwrap_or(&self.fasta_path)
                        ));
                    }
                }

                #[cfg(target_arch = "wasm32")]
                {
                    ui.vertical(|ui| {
                        ui.horizontal(|ui| {
                            if ui.button("Browse...").clicked()
                            {
                                self.open_file_picker(FilePickerType::Fasta);
                            }
                            ui.label("or drag & drop file or enter URL");
                        });

                        // Show selected file name
                        if !self.fasta_file_name.is_empty()
                        {
                            ui.label(format!("📄 File: {}", self.fasta_file_name));
                        }

                        ui.horizontal(|ui| {
                            ui.label("URL:");
                            ui.text_edit_singleline(&mut self.fasta_path);
                            if ui.button("Load").clicked()
                            {
                                self.load_fasta();
                            }
                        });
                    });
                }

                ui.separator();

                #[cfg(not(target_arch = "wasm32"))]
                {
                    if ui.button("Open GFF/GTF...").clicked()
                    {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("GFF/GTF", &["gff", "gff3", "gtf"])
                            .add_filter("GFF/GTF (gzipped)", &["gff.gz", "gff3.gz", "gtf.gz", "gz"])
                            .add_filter("All files", &["*"])
                            .pick_file()
                        {
                            self.gff_path = path.display().to_string();
                            self.load_gff();
                        }
                    }

                    if !self.gff_path.is_empty()
                    {
                        ui.label(format!(
                            "📄 {}",
                            std::path::Path::new(&self.gff_path)
                                .file_name()
                                .and_then(|n| n.to_str())
                                .unwrap_or(&self.gff_path)
                        ));
                    }
                }

                #[cfg(target_arch = "wasm32")]
                {
                    ui.vertical(|ui| {
                        ui.horizontal(|ui| {
                            if ui.button("Browse...").clicked()
                            {
                                self.open_file_picker(FilePickerType::Gff);
                            }
                            ui.label("or drag & drop file or enter URL");
                        });

                        // Show selected file name
                        if !self.gff_file_name.is_empty()
                        {
                            ui.label(format!("📄 File: {}", self.gff_file_name));
                        }

                        ui.horizontal(|ui| {
                            ui.label("URL:");
                            ui.text_edit_singleline(&mut self.gff_path);
                            if ui.button("Load").clicked()
                            {
                                self.load_gff();
                            }
                        });
                    });
                }

                ui.separator();

                // BAM file loading
                #[cfg(not(target_arch = "wasm32"))]
                {
                    if ui.button("Open BAM...").clicked()
                    {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("BAM", &["bam"])
                            .add_filter("All files", &["*"])
                            .pick_file()
                        {
                            self.bam_path = path.display().to_string();
                            self.load_bam();
                        }
                    }

                    if !self.bam_path.is_empty()
                    {
                        ui.label(format!(
                            "📊 {}",
                            std::path::Path::new(&self.bam_path)
                                .file_name()
                                .and_then(|n| n.to_str())
                                .unwrap_or(&self.bam_path)
                        ));
                    }
                }

                #[cfg(target_arch = "wasm32")]
                {
                    ui.vertical(|ui| {
                        ui.horizontal(|ui| {
                            if ui.button("Browse...").clicked()
                            {
                                self.open_file_picker(FilePickerType::Bam);
                            }
                            ui.label("BAM or drag & drop file or enter URL");
                        });

                        ui.horizontal(|ui| {
                            ui.add(
                                egui::TextEdit::singleline(&mut self.bam_path)
                                    .hint_text("https://example.com/file.bam")
                                    .desired_width(400.0),
                            );
                            if ui.button("Load").clicked()
                            {
                                self.load_bam();
                            }
                        });

                        if !self.bam_file_name.is_empty()
                        {
                            ui.label(format!("📊 {}", self.bam_file_name));
                        }
                    });
                }
            });
        });

        // Search panel (below top panel)
        egui::TopBottomPanel::top("search_panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Go to:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.position_search)
                        .hint_text("chr1:100000 or 100000")
                        .desired_width(150.0),
                );
                if ui.button("🔍 Jump").clicked() || ui.input(|i| i.key_pressed(egui::Key::Enter))
                {
                    self.search_position();
                }

                ui.separator();

                ui.label("Find feature:");
                ui.add(
                    egui::TextEdit::singleline(&mut self.feature_search)
                        .hint_text("gene name, ID...")
                        .desired_width(150.0),
                );
                if ui.button("🔍 Search").clicked()
                {
                    self.search_features();
                }
                if !self.search_results.is_empty()
                {
                    ui.label(format!("({} results)", self.search_results.len()));
                    if ui.button("Show Results").clicked()
                    {
                        self.show_search_results = !self.show_search_results;
                    }
                }

                ui.separator();

                ui.checkbox(&mut self.show_amino_acids, "Show amino acids (6 frames)");

                // BAM display toggles
                if self.alignments.is_some()
                {
                    ui.separator();
                    ui.label("BAM Tracks:");
                    ui.checkbox(&mut self.show_coverage, "Show coverage");
                    ui.checkbox(&mut self.show_alignments, "Show reads");
                    ui.checkbox(&mut self.show_variants, "Show variants");
                }
            });
        });

        egui::TopBottomPanel::bottom("bottom_panel").show(ctx, |ui| {
            ui.vertical(|ui| {
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

                // Show loading progress if active
                if let Some(ref progress) = self.loading_progress
                {
                    ui.separator();
                    ui.horizontal(|ui| {
                        ui.label(&progress.description);
                        if let Some(total) = progress.total_bytes
                        {
                            if total > 0
                            {
                                let progress_fraction = progress.bytes_loaded as f32 / total as f32;
                                ui.add(
                                    egui::ProgressBar::new(progress_fraction)
                                        .show_percentage()
                                        .text(format!(
                                            "{} / {}",
                                            format_bytes(progress.bytes_loaded),
                                            format_bytes(total)
                                        )),
                                );
                            }
                            else
                            {
                                ui.spinner();
                            }
                        }
                        else
                        {
                            ui.spinner();
                            if progress.bytes_loaded > 0
                            {
                                ui.label(format!("{} loaded", format_bytes(progress.bytes_loaded)));
                            }
                        }
                    });
                }
            });
        });

        egui::SidePanel::left("side_panel")
            .min_width(250.0)
            .show(ctx, |ui| {
                ui.heading("Chromosomes");
                ui.separator();

                if self.genome.is_some()
                {
                    // Search box
                    ui.horizontal(|ui| {
                        ui.label("🔍");
                        ui.text_edit_singleline(&mut self.chromosome_search);
                        if ui.button("✖").clicked()
                        {
                            self.chromosome_search.clear();
                        }
                    });

                    ui.separator();

                    // Sort options
                    ui.horizontal(|ui| {
                        ui.label("Sort:");
                        if ui
                            .selectable_label(
                                self.chromosome_sort == ChromosomeSort::Natural,
                                "Natural",
                            )
                            .clicked()
                        {
                            self.chromosome_sort = ChromosomeSort::Natural;
                        }
                        if ui
                            .selectable_label(
                                self.chromosome_sort == ChromosomeSort::Alphabetical,
                                "A-Z",
                            )
                            .clicked()
                        {
                            self.chromosome_sort = ChromosomeSort::Alphabetical;
                        }
                        if ui
                            .selectable_label(self.chromosome_sort == ChromosomeSort::Size, "Size")
                            .clicked()
                        {
                            self.chromosome_sort = ChromosomeSort::Size;
                        }
                    });

                    ui.separator();

                    // Chromosome list
                    let chromosomes = self.get_sorted_chromosomes();

                    if chromosomes.is_empty()
                    {
                        ui.label("No chromosomes found");
                    }
                    else
                    {
                        ui.label(format!("{} chromosome(s)", chromosomes.len()));
                        egui::ScrollArea::vertical().show(ui, |ui| {
                            for (chr_name, chr_length) in chromosomes
                            {
                                let is_selected =
                                    self.selected_chromosome.as_ref() == Some(&chr_name);
                                if ui
                                    .selectable_label(
                                        is_selected,
                                        format!("{} ({} bp)", chr_name, chr_length),
                                    )
                                    .clicked()
                                {
                                    self.selected_chromosome = Some(chr_name.clone());
                                    self.viewport = viewport::Viewport::new(0, chr_length);
                                }
                            }
                        });
                    }
                }
                else
                {
                    ui.label("No genome loaded");
                }
            });

        // Search results window
        let mut keep_results_open = self.show_search_results;
        let mut jump_to: Option<(String, usize, usize)> = None;
        let mut close_clicked = false;

        if self.show_search_results
        {
            egui::Window::new("Search Results")
                .open(&mut keep_results_open)
                .default_width(400.0)
                .show(ctx, |ui| {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        for (chr, name, start, end) in &self.search_results
                        {
                            ui.horizontal(|ui| {
                                if ui.button("Jump").clicked()
                                {
                                    jump_to = Some((chr.clone(), *start, *end));
                                }
                                ui.label(format!(
                                    "{} ({}:{}-{}, {} bp)",
                                    name,
                                    chr,
                                    start,
                                    end,
                                    end - start
                                ));
                            });
                        }
                    });

                    ui.separator();
                    if ui.button("Close").clicked()
                    {
                        close_clicked = true;
                    }
                });
        }

        if close_clicked
        {
            keep_results_open = false;
        }

        self.show_search_results = keep_results_open;

        if let Some((chr, start, end)) = jump_to
        {
            self.jump_to_feature(&chr, start, end);
        }

        egui::CentralPanel::default().show(ctx, |ui| {
            // Show drag-and-drop hint when hovering with files
            #[cfg(target_arch = "wasm32")]
            {
                ctx.input(|i| {
                    if i.raw.hovered_files.len() > 0
                    {
                        let painter = ctx.layer_painter(egui::LayerId::new(
                            egui::Order::Foreground,
                            egui::Id::new("file_drop_overlay"),
                        ));

                        let screen_rect = ctx.screen_rect();
                        painter.rect_filled(screen_rect, 0.0, egui::Color32::from_black_alpha(192));

                        painter.text(
                            screen_rect.center(),
                            egui::Align2::CENTER_CENTER,
                            "Drop FASTA or GFF/GTF file here",
                            egui::FontId::proportional(32.0),
                            egui::Color32::WHITE,
                        );
                    }
                });
            }

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
                // Horizontal scroll (two-finger swipe left/right) for pan
                let scroll_delta_x = ui.input(|i| i.smooth_scroll_delta.x);

                // Vertical scroll (two-finger swipe up/down) for zoom
                let scroll_delta_y = ui.input(|i| i.smooth_scroll_delta.y);

                if scroll_delta_y > scroll_delta_x
                {
                    if scroll_delta_y != 0.0
                    {
                        let mouse_x = hover_pos.x - response.rect.left();
                        let relative_pos = mouse_x / available_size.x;
                        let zoom_factor = if scroll_delta_y > 0.0 { 0.9 } else { 1.1 };
                        self.viewport.zoom(zoom_factor, relative_pos);
                    }
                }
                else
                {
                    if scroll_delta_x != 0.0
                    {
                        let pixels_per_base =
                            available_size.x / (self.viewport.end - self.viewport.start) as f32;
                        let bases_moved = (scroll_delta_x / pixels_per_base) as i64;
                        self.viewport.pan(bases_moved);
                    }
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
                        self.show_amino_acids,
                    );

                    // Render BAM tracks (if loaded and selected chromosome matches)
                    if let Some(ref alignments) = self.alignments
                    {
                        if let Some(track) = alignments.tracks.get(chr_name)
                        {
                            // Start BAM tracks below genome tracks
                            // Estimate y_offset based on visible genome tracks
                            let mut bam_y_offset = response.rect.top() + 150.0;

                            // Adjust based on what's visible
                            if self.viewport.width() < 1000
                            {
                                bam_y_offset += 50.0; // sequence track
                            }
                            if self.show_amino_acids && self.viewport.width() < 5000
                            {
                                bam_y_offset += 130.0; // amino acid frames
                            }
                            if !self.features.is_empty()
                            {
                                bam_y_offset += 100.0; // feature tracks
                            }

                            if self.show_coverage
                            {
                                bam_y_offset = renderer::draw_coverage_track(
                                    &painter,
                                    response.rect,
                                    track,
                                    &self.viewport,
                                    bam_y_offset,
                                );
                            }

                            if self.show_alignments
                            {
                                bam_y_offset = renderer::draw_alignments(
                                    &painter,
                                    response.rect,
                                    track,
                                    &self.viewport,
                                    bam_y_offset,
                                );
                            }

                            if self.show_variants
                            {
                                renderer::draw_variant_summary(
                                    &painter,
                                    response.rect,
                                    track,
                                    &self.viewport,
                                    bam_y_offset,
                                );
                            }
                        }
                    }
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
