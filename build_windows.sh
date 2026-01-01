#!/bin/bash
set -e

echo "Building UGV for Windows (x86_64)..."

# Check if Windows target is installed
if ! rustup target list --installed | grep -q "x86_64-pc-windows-gnu"; then
    echo "Installing Windows target..."
    rustup target add x86_64-pc-windows-gnu
fi

# Check if mingw-w64 is installed
if ! command -v x86_64-w64-mingw32-gcc &> /dev/null; then
    echo ""
    echo "ERROR: mingw-w64 is not installed!"
    echo "Please install it first:"
    echo "  Ubuntu/Debian: sudo apt install mingw-w64"
    echo "  Fedora:        sudo dnf install mingw64-gcc"
    echo "  Arch:          sudo pacman -S mingw-w64-gcc"
    exit 1
fi

# Build for Windows
echo "Building release binary..."
cargo build --release --target x86_64-pc-windows-gnu

echo ""
echo "Build complete!"
echo "Windows executable: target/x86_64-pc-windows-gnu/release/ugv.exe"
