# UGV - Ultra-Fast Genome Viewer

A high-performance genome viewer built with Rust and egui, designed for interactive exploration of genomic data.

## Features

### Genome Visualization
- **FASTA Support**: Load genome sequences in FASTA format (plain or gzipped)
- **GFF/GTF Annotations**: Display gene features from GFF3/GTF files (plain or gzipped)
- **Multi-track Display**:
  - Position ruler with adaptive scaling
  - GC content plot
  - DNA sequence view (when zoomed in)
  - **Amino acid translation** in all 6 reading frames (3 forward + 3 reverse)
  - Gene feature tracks with color-coding by type

### Interactive Navigation
- **Pan**:
  - Click and drag to move along the chromosome
  - Two-finger horizontal swipe (touchpad/trackpad) to scroll sideways
- **Zoom**:
  - Mouse scroll wheel to zoom in/out (focus-aware zooming)
  - Two-finger vertical swipe (touchpad/trackpad) to zoom
- **Multi-chromosome Support**: Browse all chromosomes in the genome

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

## Installation

### Prerequisites
- Rust 1.70 or later

### Build from Source
```bash
git clone <repository-url>
cd ugv
cargo build --release
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

3. **Navigate**:
   - Select a chromosome from the left panel
   - Use the search box to filter chromosomes
   - Click sort buttons to change chromosome order
   - Drag to pan, scroll to zoom

4. **Search and Jump**:
   - **Go to position**: Enter `chr1:100000` in the "Go to:" field and press Enter
   - **Find features**: Enter a gene name in "Find feature:" and click Search
   - Browse results and click "Jump" to navigate to any feature

5. **View Amino Acid Translation**:
   - Enable the "Show amino acids (6 frames)" checkbox in the search panel
   - Zoom in to view level (< 5000 bases) to see the amino acid translations
   - All 6 reading frames are displayed (3 forward + 3 reverse)
   - Amino acids are color-coded by biochemical properties

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

## Navigation Controls
- **Mouse wheel / Two-finger vertical swipe**: Zoom in/out (focus-aware)
- **Click + drag**: Pan view horizontally
- **Two-finger horizontal swipe**: Scroll genome sideways
- **Enter** (in position search): Jump to position

## Architecture

### Modules
- **fasta.rs**: FASTA genome parser with GC content calculation
- **gff.rs**: GFF3/GTF annotation parser
- **interval_tree.rs**: Efficient feature range queries
- **viewport.rs**: View management (pan, zoom, coordinate mapping)
- **translation.rs**: DNA to protein translation (standard genetic code, 6 frames)
- **renderer.rs**: Multi-track genome visualization with amino acid display

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
1. Drag a FASTA or GFF/GTF file from your computer onto the browser window
2. The file will be automatically detected by extension and loaded
3. Supports both plain and gzipped files (`.gz`)
4. A visual overlay appears when hovering with files

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
