use anyhow::{Context, Result};
use std::io::Read;

/// Load file content from either a local path or HTTP URL
pub fn load_file(path: &str) -> Result<Vec<u8>> {
    #[cfg(not(target_arch = "wasm32"))]
    {
        load_file_native(path)
    }

    #[cfg(target_arch = "wasm32")]
    {
        Err(anyhow::anyhow!(
            "Synchronous file loading not supported in WASM. Use load_file_async instead."
        ))
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn load_file_native(path: &str) -> Result<Vec<u8>> {
    use std::fs::File;
    use std::io::BufReader;

    // Check if it's a URL
    if path.starts_with("http://") || path.starts_with("https://") {
        return Err(anyhow::anyhow!(
            "HTTP URLs not supported in native build. Please download the file first."
        ));
    }

    let file = File::open(path).with_context(|| format!("Failed to open file: {}", path))?;
    let mut reader = BufReader::new(file);
    let mut buffer = Vec::new();
    reader
        .read_to_end(&mut buffer)
        .with_context(|| format!("Failed to read file: {}", path))?;
    Ok(buffer)
}

#[cfg(target_arch = "wasm32")]
pub async fn load_file_async(path: &str) -> Result<Vec<u8>> {
    use gloo_net::http::Request;

    if !path.starts_with("http://") && !path.starts_with("https://") {
        return Err(anyhow::anyhow!(
            "In WASM, only HTTP/HTTPS URLs are supported. Got: {}",
            path
        ));
    }

    let response = Request::get(path)
        .send()
        .await
        .map_err(|e| anyhow::anyhow!("Failed to fetch {}: {}", path, e))?;

    if !response.ok() {
        return Err(anyhow::anyhow!(
            "HTTP error {}: {}",
            response.status(),
            response.status_text()
        ));
    }

    let bytes = response
        .binary()
        .await
        .map_err(|e| anyhow::anyhow!("Failed to read response body: {}", e))?;

    Ok(bytes)
}
