# UGV - Ultra-Fast Genome Viewer

A high-performance genome viewer built with Rust and egui, designed for interactive exploration of genomic data.

## Features

### Genome Visualization
- **FASTA Support**: Load genome sequences in FASTA format (plain or gzipped)
- **GFF/GTF Annotations**: Display gene features from GFF3/GTF files (plain or gzipped)
- **BAM Sequencing Data**: Display aligned sequencing reads from BAM files
- **Custom TSV Tracks**: Load and visualize quantitative genomic data from TSV files
- **Multi-track Display**:
  - Position ruler with adaptive scaling
  - GC content plot
  - DNA sequence view (when zoomed in)
  - **Amino acid translation** in all 6 reading frames (3 forward + 3 reverse)
  - Gene feature tracks with color-coding by type
  - **Sequencing tracks**:
    - Coverage histogram showing read depth
    - Aligned reads with stacking to avoid overlaps
    - Variant summary highlighting SNPs, insertions, and deletions
  - **Custom data tracks**:
    - Line graph visualization of quantitative data
    - Auto-scaling to fit data range
    - Zero line indicator
    - Position-based plotting

### Interactive Navigation
- **Pan**:
  - Click and drag to move along the chromosome
  - Two-finger horizontal swipe (touchpad/trackpad) to scroll sideways
- **Zoom**:
  - Mouse scroll wheel to zoom in/out (focus-aware zooming)
  - Two-finger vertical swipe (touchpad/trackpad) to zoom
- **Multi-chromosome Support**: Browse all chromosomes in the genome

### Session Management
- **Save Sessions**: Save your current workspace including:
  - All loaded file paths (FASTA, GFF, BAM)
  - Current chromosome and viewport position
  - UI settings (chromosome panel visibility, amino acid display, BAM track toggles)
  - Chromosome sorting preference
- **Load Sessions**: Restore your workspace from a saved session file
- **Session Format**: JSON format (.ugv or .json extension)
- **Cross-platform**: Works in both native and WebAssembly builds
  - Native: File dialog for save/load
  - WASM: Downloads session file / loads from file picker or drag & drop

### Chromosome Management
- **Search**: Filter chromosomes by name
- **Sorting Options**:
  - **Natural** (default): chr1, chr2, ..., chr10, ..., chrX, chrY, chrM
  - **Alphabetical**: Sort by name A-Z
  - **Size**: Sort by chromosome length (largest first)

### Navigation and Search
- **Position Search**: Jump to specific genomic positions
  - Format: `chr1:100000` or just `100000` (uses current chromosome)
  - Supports comma/underscore separators: `chr1:100,000` or `chr1:100_000`
  - Press Enter or click "Jump" to navigate
- **Feature Search**: Find genes and annotations by name
  - Search by gene name, ID, or any GFF3 attribute
  - Case-insensitive partial matching
  - View all results in a popup window
  - Click "Jump" on any result to navigate to that feature

### Performance Features
- Interval tree for efficient feature queries (O(log n))
- On-the-fly gzip decompression
- GPU-accelerated rendering via egui
- Responsive viewport clipping

### Feature Color Coding
- **Gene**: Steel Blue
- **mRNA/Transcript**: Light Blue
- **Exon**: Forest Green
- **CDS**: Orange
- **UTR**: Yellow
- **Intron**: Gray

### Amino Acid Color Coding
When amino acid frames are enabled, amino acids are color-coded by type:
- **Hydrophobic (A, V, I, L, M)**: Green
- **Aromatic (F, W, Y)**: Purple
- **Polar (S, T, N, Q)**: Teal
- **Positively charged (K, R, H)**: Blue
- **Negatively charged (D, E)**: Red
- **Cysteine (C)**: Yellow
- **Glycine (G)**: Gray
- **Proline (P)**: Orange
- **Stop codon (*)**: Red

#### Reading Frames
- **Forward frames** (+1, +2, +3): Light blue background
- **Reverse frames** (-1, -2, -3): Light orange background

### Sequencing Track Color Coding

#### Aligned Reads
- **Forward strand**: Blue
- **Reverse strand**: Red
- **Opacity**: Based on mapping quality (higher quality = more opaque)

#### Variant Summary
- **SNPs (Single Nucleotide Polymorphisms)**: Red
- **Insertions**: Purple
- **Deletions**: Black

## Installation

### Prerequisites
- Rust 1.70 or later

### Build from Source

#### Native Build (Linux/macOS/Windows)
```bash
git clone <repository-url>
cd ugv
cargo build --release
```

#### Cross-Compile for Windows (from Linux)
```bash
./build_windows.sh
```

Requirements:
- `mingw-w64` toolchain
  - Ubuntu/Debian: `sudo apt install mingw-w64`
  - Fedora: `sudo dnf install mingw64-gcc`
  - Arch: `sudo pacman -S mingw-w64-gcc`

Output: `target/x86_64-pc-windows-gnu/release/ugv.exe`

#### WebAssembly Build
```bash
./build_wasm.sh
```

This will:
- Build the WASM binary
- Generate JavaScript bindings
- Create files ready for web deployment

To test locally:
```bash
python3 -m http.server 8080
# Open http://localhost:8080 in your browser
```

### Run
```bash
cargo run --release
```

## Usage

1. **Load Genome**: Click "Open FASTA..." and select your genome file
   - Supports: `.fasta`, `.fa`, `.fna`, `.ffn`, `.faa`, `.frn`
   - Gzipped: `.fasta.gz`, `.fa.gz`, etc.

2. **Load Annotations** (optional): Click "Open GFF/GTF..." and select your annotation file
   - Supports: `.gff`, `.gff3`, `.gtf`
   - Gzipped: `.gff.gz`, `.gff3.gz`, `.gtf.gz`

3. **Load Sequencing Data** (optional): Click "Open BAM..." and select your alignment file
   - Supports: `.bam` (Binary Alignment/Map format)
   - Works with both native and WebAssembly builds
   - Enable/disable individual tracks using checkboxes:
     - "Show coverage": Read depth histogram
     - "Show reads": Individual aligned reads with stacking
     - "Show variants": SNPs and indels detected from CIGAR strings

4. **Load Custom Data Track** (optional): Click "Open TSV..." and select your data file
   - Supports: `.tsv`, `.txt` (Tab-separated values)
   - Format: chromosome TAB position TAB signal-value
   - Displays quantitative genomic data as a line graph
   - Auto-scales to data range in current view
   - Enable/disable using "Show custom TSV track" checkbox

5. **Navigate**:
   - Select a chromosome from the left panel
   - Use the search box to filter chromosomes
   - Click sort buttons to change chromosome order
   - Drag to pan, scroll to zoom

6. **Search and Jump**:
   - **Go to position**: Enter `chr1:100000` in the "Go to:" field and press Enter
   - **Find features**: Enter a gene name in "Find feature:" and click Search
   - Browse results and click "Jump" to navigate to any feature

7. **View Amino Acid Translation**:
   - Enable the "Show amino acids (6 frames)" checkbox in the search panel
   - Zoom in to view level (< 5000 bases) to see the amino acid translations
   - All 6 reading frames are displayed (3 forward + 3 reverse)
   - Amino acids are color-coded by biochemical properties

8. **Save and Load Sessions**:
   - Click **💾 Save** in the Session section to save your current workspace
   - Saves file paths, viewport position, and all UI settings to a .ugv file
   - Click **📂 Load** to restore a previously saved session
   - Sessions can be shared or used to quickly return to specific analyses
   - Supports drag & drop for session files (WASM version)

## File Format Support

### FASTA Files
Standard FASTA format for genome sequences:
```
>chr1
ATCGATCGATCG...
>chr2
GCTAGCTAGCTA...
```

### GFF3/GTF Files
Standard GFF3 or GTF format for gene annotations:
```
chr1    source    gene    1000    5000    .    +    .    ID=gene1;Name=MyGene
chr1    source    exon    1000    1500    .    +    .    Parent=gene1
```

### BAM Files
Binary Alignment/Map (BAM) format for aligned sequencing reads:
- Stores aligned reads from sequencing experiments (DNA-seq, RNA-seq, etc.)
- Binary compressed format for efficient storage
- Parses CIGAR strings to detect variants (SNPs, insertions, deletions)
- Computes coverage histograms with 100bp binning for performance
- Displays up to 50 rows of stacked reads to avoid overlaps
- Color-coded by strand (forward/reverse) and mapping quality

**Supported Features:**
- Read alignment visualization with CIGAR parsing
- Coverage depth calculation across the genome
- Variant detection from CIGAR operations (Match, Insertion, Deletion, Skip)
- Strand-specific display (forward/reverse)
- Mapping quality visualization (alpha blending)

**Performance Notes:**
- Coverage is pre-computed in 100bp bins for fast rendering
- At >1000 visible reads, shows warning message and coverage only
- Maximum 50 rows of stacked reads displayed
- Works in both native and WebAssembly builds

### TSV Custom Track Files
Tab-separated values format for quantitative genomic data:

**Format:**
```tsv
chr1	100	1.5
chr1	200	2.3
chr1	300	-0.5
chr2	150	3.2
```

**Field Specifications:**
- **Field 1** (chromosome): Chromosome name (must match FASTA chromosome names)
- **Field 2** (position): Genomic position (0-based or 1-based, integer)
- **Field 3** (signal): Quantitative value (float64)

**Features:**
- Simple tab-delimited format
- File extensions: `.tsv` or `.txt`
- Comments supported (lines starting with `#`)
- Empty lines ignored
- Flexible: works with any quantitative genomic data
- Visualized as line graph with auto-scaling
- Examples: expression levels, conservation scores, methylation ratios, etc.

**Example Use Cases:**
- Gene expression levels across genome
- Conservation scores (PhyloP, PhastCons)
- ChIP-seq signal tracks
- Methylation percentages
- Custom scoring metrics

### Session Files
Session files store your workspace state in JSON format:
```json
{
  "fasta_path": "/path/to/genome.fasta",
  "gff_path": "/path/to/annotations.gff",
  "bam_path": "/path/to/alignments.bam",
  "selected_chromosome": "chr1",
  "viewport_start": 100000,
  "viewport_end": 200000,
  "show_chromosome_panel": true,
  "show_amino_acids": false,
  "show_coverage": true,
  "show_alignments": true,
  "show_variants": false,
  "chromosome_sort": "Natural"
}
```

**Features:**
- Human-readable JSON format
- File extensions: `.ugv` (recommended) or `.json`
- Stores absolute or relative file paths
- Automatically loads all referenced files when session is restored
- UI settings preserved (panel visibility, display options, sort order)
- Can be edited manually if needed

## Navigation Controls
- **Mouse wheel / Two-finger vertical swipe**: Zoom in/out (focus-aware)
- **Click + drag**: Pan view horizontally
- **Two-finger horizontal swipe**: Scroll genome sideways
- **Enter** (in position search): Jump to position

## Architecture

### Modules
- **fasta.rs**: FASTA genome parser with GC content calculation
- **gff.rs**: GFF3/GTF annotation parser
- **bam.rs**: BAM/SAM parser with CIGAR operations, variant extraction, and coverage calculation
- **tsv.rs**: TSV custom track parser for quantitative genomic data
- **session.rs**: Session management for saving and restoring workspace state
- **interval_tree.rs**: Efficient feature range queries
- **viewport.rs**: View management (pan, zoom, coordinate mapping)
- **translation.rs**: DNA to protein translation (standard genetic code, 6 frames)
- **renderer.rs**: Multi-track genome visualization with amino acid display, sequencing tracks, and custom data tracks

### Performance Optimizations
- Binary search for interval queries
- Viewport-based rendering (only visible features)
- Adaptive ruler tick intervals
- Efficient GC content windowing

## WebAssembly Support

The viewer can be compiled to WebAssembly for browser deployment.

### Quick Build

Use the provided build script:

```bash
chmod +x build_wasm.sh
./build_wasm.sh
```

This will:
1. Build the project for wasm32
2. Install `wasm-bindgen-cli` if needed
3. Generate JavaScript bindings

### Run Locally

Start a local web server:

```bash
python3 -m http.server 8080
```

Then open http://localhost:8080 in your browser.

### Manual Build

If you prefer to build manually:

```bash
# Install wasm toolchain
rustup target add wasm32-unknown-unknown

# Build for web
cargo build --release --target wasm32-unknown-unknown

# Install wasm-bindgen-cli
cargo install wasm-bindgen-cli

# Generate bindings
wasm-bindgen target/wasm32-unknown-unknown/release/ugv.wasm \
    --out-dir . \
    --target web \
    --no-typescript
```

### WASM Limitations

Due to browser security restrictions, the WebAssembly version:
- Uses text input fields instead of native file dialogs
- Requires files to be accessible via HTTP/HTTPS URLs
- Automatically fetches and loads files from provided URLs
- Supports gzipped files (`.gz`) with automatic decompression

### Using the WASM Version

The WebAssembly version supports three methods for loading files:

#### Method 1: File Dialog (Most User-Friendly)
1. Click the "Browse..." button for FASTA or GFF/GTF files
2. Use your browser's native file picker to select a local file
3. The file will be automatically loaded and parsed
4. Supports both plain and gzipped files (`.gz`)

#### Method 2: Drag and Drop
1. Drag a FASTA, GFF/GTF, BAM, or session file from your computer onto the browser window
2. The file will be automatically detected by extension and loaded
3. Supports both plain and gzipped files (`.gz`)
4. A visual overlay appears when hovering with files
5. BAM files (`.bam`) and session files (`.ugv`, `.json`) are automatically recognized and parsed
6. Session files restore your complete workspace including all file paths and settings

#### Method 3: HTTP/HTTPS URL
1. Enter the full URL to your genome file in the URL field
   - Example: `https://ftp.ensembl.org/pub/release-115/fasta/bos_taurus/dna_index/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa.gz`
2. Click "Load" and wait for the file to download and parse
3. Useful for loading files from public genome databases

**Loading Progress:**
- Status bar shows current operation ("Loading FASTA from...")
- Progress bar displays loading/parsing progress with file size information
- Visual spinner indicates ongoing operations
- Shows download progress for HTTP URLs (when content-length is available)
- Shows file size and parsing status for local files
- Large files may take time to download and parse
- Gzipped files are automatically decompressed
- File name is displayed once loaded

## License

[Add your license here]

## Contributing

Contributions welcome! Please feel free to submit issues or pull requests.
