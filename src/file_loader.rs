use anyhow::Result;

/// Load file content from either a local path or HTTP URL
#[allow(unused_variables)]
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
    use anyhow::Context;
    use std::fs::File;
    use std::io::{BufReader, Read};

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
    load_file_with_progress(path, |_, _| {}).await
}

#[cfg(target_arch = "wasm32")]
pub async fn load_file_with_progress<F>(path: &str, mut progress_callback: F) -> Result<Vec<u8>>
where
    F: FnMut(usize, Option<usize>),
{
    use wasm_bindgen::JsCast;
    use wasm_bindgen_futures::JsFuture;
    use web_sys::Response;

    if !path.starts_with("http://") && !path.starts_with("https://") {
        return Err(anyhow::anyhow!(
            "In WASM, only HTTP/HTTPS URLs are supported. Got: {}",
            path
        ));
    }

    let window = web_sys::window().ok_or_else(|| anyhow::anyhow!("No window object"))?;

    let response_promise = window.fetch_with_str(path);
    let response_value = JsFuture::from(response_promise)
        .await
        .map_err(|e| anyhow::anyhow!("Fetch failed: {:?}", e))?;

    let response: Response = response_value
        .dyn_into()
        .map_err(|_| anyhow::anyhow!("Response cast failed"))?;

    if !response.ok() {
        return Err(anyhow::anyhow!(
            "HTTP error {}: {}",
            response.status(),
            response.status_text()
        ));
    }

    // Get content length if available
    let content_length = response
        .headers()
        .get("content-length")
        .ok()
        .flatten()
        .and_then(|s| s.parse::<usize>().ok());

    // Read the body
    let array_buffer_promise = response
        .array_buffer()
        .map_err(|_| anyhow::anyhow!("Failed to get array buffer"))?;

    let array_buffer_value = JsFuture::from(array_buffer_promise)
        .await
        .map_err(|e| anyhow::anyhow!("Failed to read array buffer: {:?}", e))?;

    let uint8_array = js_sys::Uint8Array::new(&array_buffer_value);
    let bytes = uint8_array.to_vec();

    // Report final progress
    progress_callback(bytes.len(), content_length);

    Ok(bytes)
}
