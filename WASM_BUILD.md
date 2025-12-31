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
web-sys = "0.3"
log = "0.4"
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

In WASM builds, users provide file URLs instead of selecting local files:

- Genome files must be served via HTTP/HTTPS
- Files can be hosted on the same origin or via CORS-enabled servers
- Example: `https://example.com/genomes/chr1.fasta.gz`

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

## Future Enhancements

Potential WASM-specific improvements:

1. **Drag & Drop**: Support dragging genome files onto the browser
2. **File API**: Use browser File API for local file access
3. **IndexedDB**: Cache loaded genomes in browser storage
4. **SharedArrayBuffer**: Enable for multi-threading (requires secure context)
5. **Streaming**: Progressive loading of large genome files

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
