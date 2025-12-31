#!/bin/bash
set -e

echo "Building UGV for WebAssembly..."

# Build for wasm32
cargo build --release --target wasm32-unknown-unknown

# Install wasm-bindgen-cli if not already installed
if ! command -v wasm-bindgen &> /dev/null; then
    echo "Installing wasm-bindgen-cli..."
    cargo install wasm-bindgen-cli
fi

# Generate JavaScript bindings
echo "Generating JavaScript bindings..."
wasm-bindgen target/wasm32-unknown-unknown/release/ugv.wasm \
    --out-dir . \
    --target web \
    --no-typescript

echo "Build complete!"
echo "To run locally, start a web server in this directory:"
echo "  python3 -m http.server 8080"
echo "Then open http://localhost:8080 in your browser"
