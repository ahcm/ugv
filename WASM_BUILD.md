# WebAssembly Build Guide

This document explains how the UGV genome viewer has been configured to support both native and WebAssembly builds.

## Changes Made for WASM Support

### 1. Conditional Compilation

The codebase uses Rust's conditional compilation features to support both platforms:

**Main Function** (src/main.rs):
- Native builds use `eframe::run_native()` with `NativeOptions`
- WASM builds use `eframe::WebRunner` with `WebOptions`

**File Loading UI**:
- Native: Uses `rfd` file dialogs
- WASM: Uses text input fields for file paths/URLs

### 2. Dependencies

**Cargo.toml** has been split into platform-specific dependencies:

```toml
[dependencies]
# Shared dependencies
eframe = "0.29"
egui = "0.29"
# ... other shared deps

[target.'cfg(not(target_arch = "wasm32"))'.dependencies]
rfd = "0.15"  # Native file dialogs only

[target.'cfg(target_arch = "wasm32")'.dependencies]
wasm-bindgen = "0.2"
wasm-bindgen-futures = "0.4"
js-sys = "0.3"            # JavaScript interop (Uint8Array, etc.)
web-sys = { version = "0.3", features = [
    "Request", "RequestInit", "RequestMode", "Response", "Window",
    "File", "FileList", "FileReader", "Blob", "Headers",
    "Document", "Element", "HtmlInputElement", "EventTarget", "Event"
] }
log = "0.4"
gloo-net = "0.6"          # HTTP requests
gloo-file = "0.3"         # File reading
```

### 3. Build Process

The `build_wasm.sh` script automates the WASM build:

1. Compiles to `wasm32-unknown-unknown` target
2. Uses `wasm-bindgen` to generate JavaScript bindings
3. Produces `ugv.js` and `ugv_bg.wasm` files

### 4. Web Deployment

**index.html** provides the HTML shell for the WASM app:
- Creates the canvas element with id "the_canvas_id"
- Loads and initializes the WASM module
- Shows loading spinner during initialization

## Building for WebAssembly

### Prerequisites

```bash
# Install Rust wasm32 target
rustup target add wasm32-unknown-unknown

# Install wasm-bindgen-cli (done automatically by build_wasm.sh)
cargo install wasm-bindgen-cli
```

### Build Commands

**Quick build:**
```bash
./build_wasm.sh
```

**Manual build:**
```bash
cargo build --release --target wasm32-unknown-unknown
wasm-bindgen target/wasm32-unknown-unknown/release/ugv.wasm \
    --out-dir . \
    --target web \
    --no-typescript
```

### Running Locally

```bash
# Start a local server
python3 -m http.server 8080

# Open browser to http://localhost:8080
```

## Platform Differences

### Native Build
- **File Loading**: Native OS file picker dialogs
- **File Access**: Direct filesystem access
- **Performance**: Full native performance
- **File Formats**: Can read local .gz files directly

### WASM Build
- **File Loading**: Text input fields for URLs/paths
- **File Access**: HTTP/HTTPS URLs or data URLs only
- **Performance**: Near-native with WASM optimizations
- **File Formats**: Supports .gz decompression in browser

## Technical Details

### WASM Initialization

The WASM version uses async initialization:

```rust
wasm_bindgen_futures::spawn_local(async {
    let canvas = document.get_element_by_id("the_canvas_id")...;
    eframe::WebRunner::new()
        .start(canvas, web_options, Box::new(|_cc| Ok(Box::new(GenomeViewer::new()))))
        .await
        .expect("failed to start eframe");
});
```

### File Loading in WASM

The WASM build supports multiple file loading methods:

#### HTTP/HTTPS URLs (Recommended for FASTA files)
- **FASTA files MUST be loaded from URLs with index files (.fai and .gzi)**
- Index files enable efficient chromosome loading without downloading entire genomes
- Files can be hosted on the same origin or via CORS-enabled servers
- Example: `https://ftp.ensembl.org/pub/release-115/fasta/bos_taurus/dna_index/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa.gz`
- Required index files:
  - `.fai` - FASTA index (always required)
  - `.gzi` - BGZF index (required for .gz files)
- Useful for loading from public genome databases like Ensembl, NCBI

#### File Dialog (Browser Native)
- "Browse..." button opens the browser's native file picker
- **For FASTA files**: Not supported - use URLs instead (prevents browser hang on large files)
- **For GFF/GTF/BAM/TSV files**: Supported for local file upload
- Uses HTML5 File API via web-sys bindings
- File type filtering (.gff, .bam, .tsv, .gz extensions)
- Automatic file type detection and loading

#### Drag and Drop
- Users can drag local files directly onto the browser window
- **For FASTA files**: Not supported - use URLs instead
- **For GFF/GTF/BAM/TSV files**: Supported for local file upload
- Files are read using the browser's File API
- Automatic file type detection based on extension
- Visual feedback with overlay when hovering with files
- No need to upload files to a server

#### File Format Support
- FASTA files: `.fasta`, `.fa`, `.fna`, `.ffn`, `.faa`, `.frn`
- GFF/GTF files: `.gff`, `.gff3`, `.gtf`
- All formats support gzip compression (`.gz` extension)

#### Loading Progress Tracking
The WASM version provides visual feedback during file loading:

- **Progress Bar**: Shows percentage complete for local files (file size known)
- **Spinner**: Animated spinner for HTTP downloads (when content-length unavailable)
- **Status Messages**: Updates showing current operation:
  - "Downloading FASTA file..." - Fetching from URL
  - "Reading FASTA file..." - Loading local file
  - "Parsing FASTA..." - Processing file contents
- **File Size Display**: Shows bytes loaded and total file size in human-readable format (KB, MB, GB)
- **Bottom Panel**: Progress information displayed in the bottom status panel

## Deployment Options

### Static Hosting

Deploy to any static hosting service:
- GitHub Pages
- Netlify
- Vercel
- AWS S3 + CloudFront
- Any web server

Required files:
- `index.html`
- `ugv.js`
- `ugv_bg.wasm`

### CORS Considerations

If loading genome files from external sources, ensure:
- The hosting server has CORS headers configured
- `Access-Control-Allow-Origin` permits your domain

## Creating Index Files for FASTA

The web version requires indexed FASTA files for efficient loading. Here's how to create the required index files:

### Using samtools

```bash
# For uncompressed FASTA
samtools faidx genome.fa

# For BGZF-compressed FASTA (recommended for web use)
# 1. Compress with bgzip (NOT regular gzip)
bgzip genome.fa  # Creates genome.fa.gz

# 2. Create FASTA index
samtools faidx genome.fa.gz  # Creates genome.fa.gz.fai

# 3. Create BGZF index
samtools gzi genome.fa.gz    # Creates genome.fa.gz.gzi
```

### Required Files

After indexing, you need these files:
- `genome.fa.gz` - BGZF-compressed FASTA
- `genome.fa.gz.fai` - FASTA index (5-column tab-separated file)
- `genome.fa.gz.gzi` - BGZF index (binary format)

### Hosting the Files

Upload all three files to a web server with CORS enabled:
```
https://example.com/genomes/genome.fa.gz
https://example.com/genomes/genome.fa.gz.fai
https://example.com/genomes/genome.fa.gz.gzi
```

Then use the base URL in UGV:
```
https://example.com/genomes/genome.fa.gz
```

The index files will be automatically fetched.

## Implemented Features

WASM-specific features already implemented:

1. ✅ **Indexed FASTA Loading**: Efficient chromosome-level loading using .fai and .gzi index files
2. ✅ **HTTP Loading**: Load files from URLs (Ensembl, NCBI, etc.)
3. ✅ **File Dialog**: Native browser file picker for GFF/BAM/TSV files
4. ✅ **Drag & Drop**: Drag GFF/BAM/TSV files onto the browser window
5. ✅ **File API**: Uses browser File API for local file access (non-FASTA files)
6. ✅ **Gzip Support**: Automatic decompression of `.gz` files with BGZF indexing
7. ✅ **Visual Feedback**: Overlay when hovering with files
8. ✅ **Loading Progress**: Real-time progress bar showing file loading and parsing status
9. ✅ **Touch/Trackpad Gestures**: Two-finger swipe horizontal for panning, vertical for zooming
10. ✅ **Chromosome Caching**: Lazy-load and cache chromosomes on demand

## Future Enhancements

Potential improvements:

1. **IndexedDB**: Cache loaded genomes in browser storage for faster reload
2. **SharedArrayBuffer**: Enable for multi-threading (requires secure context)
3. **Streaming**: Progressive loading of large genome files
4. **Service Worker**: Offline support and background file processing
5. **Chunked Loading**: Load only visible chromosome regions for very large genomes

## Troubleshooting

### Build Errors

**Error: `wasm-bindgen` version mismatch**
```bash
cargo install wasm-bindgen-cli --force
```

**Error: Cannot find module**
- Ensure `ugv.js` and `ugv_bg.wasm` are in the same directory as `index.html`
- Check browser console for CORS errors

### Runtime Errors

**Canvas not found**
- Verify `index.html` has `<canvas id="the_canvas_id">`

**File loading fails**
- Check file URLs are accessible
- Verify CORS headers if loading from external domain
- Test with browser developer tools network tab

## References

- [eframe Web Documentation](https://docs.rs/eframe/latest/eframe/)
- [wasm-bindgen Guide](https://rustwasm.github.io/wasm-bindgen/)
- [Rust and WebAssembly Book](https://rustwasm.github.io/docs/book/)
