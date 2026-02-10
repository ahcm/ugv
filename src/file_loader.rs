use anyhow::Result;

/// Load file content from either a local path or HTTP URL
#[allow(unused_variables)]
pub fn load_file(path: &str) -> Result<Vec<u8>>
{
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
fn load_file_native(path: &str) -> Result<Vec<u8>>
{
    use anyhow::Context;
    use std::fs::File;
    use std::io::{BufReader, Read};

    // Check if it's a URL
    if path.starts_with("http://") || path.starts_with("https://")
    {
        // Use reqwest to download the file
        let response = reqwest::blocking::get(path)
            .with_context(|| format!("Failed to download from URL: {}", path))?;

        if !response.status().is_success()
        {
            return Err(anyhow::anyhow!(
                "HTTP error {}: Failed to download {}",
                response.status(),
                path
            ));
        }

        let bytes = response
            .bytes()
            .with_context(|| format!("Failed to read response body from: {}", path))?;

        return Ok(bytes.to_vec());
    }

    // Local file path
    let file = File::open(path).with_context(|| format!("Failed to open file: {}", path))?;
    let mut reader = BufReader::new(file);
    let mut buffer = Vec::new();
    reader
        .read_to_end(&mut buffer)
        .with_context(|| format!("Failed to read file: {}", path))?;
    Ok(buffer)
}

#[cfg(target_arch = "wasm32")]
pub async fn load_file_async(path: &str) -> Result<Vec<u8>>
{
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

    if !path.starts_with("http://") && !path.starts_with("https://")
    {
        return Err(anyhow::anyhow!("In WASM, only HTTP/HTTPS URLs are supported. Got: {}", path));
    }

    let window = web_sys::window().ok_or_else(|| anyhow::anyhow!("No window object"))?;

    let response_promise = window.fetch_with_str(path);
    let response_value = JsFuture::from(response_promise)
        .await
        .map_err(|e| anyhow::anyhow!("Fetch failed: {:?}", e))?;

    let response: Response = response_value
        .dyn_into()
        .map_err(|_| anyhow::anyhow!("Response cast failed"))?;

    if !response.ok()
    {
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

/// Check if index files (.fai and .gzi) exist for a given FASTA URL
/// Returns (fai_url, gzi_url) if both exist, None otherwise
#[cfg(target_arch = "wasm32")]
pub async fn check_index_files(fasta_url: &str) -> Option<(String, String)>
{
    let fai_url = format!("{}.fai", fasta_url);
    let gzi_url = format!("{}.gzi", fasta_url);

    // Check if FAI exists
    let fai_exists = check_url_exists(&fai_url).await;
    if !fai_exists
    {
        return None;
    }

    // Check if GZI exists (only needed for gzipped files)
    if fasta_url.ends_with(".gz")
    {
        let gzi_exists = check_url_exists(&gzi_url).await;
        if !gzi_exists
        {
            return None;
        }
    }

    Some((fai_url, gzi_url))
}

/// Check if a URL returns a successful response (HEAD or GET request)
#[cfg(target_arch = "wasm32")]
async fn check_url_exists(url: &str) -> bool
{
    use wasm_bindgen::JsCast;
    use wasm_bindgen_futures::JsFuture;
    use web_sys::{RequestInit, RequestMode, Response};

    let window = match web_sys::window()
    {
        Some(w) => w,
        None => return false,
    };

    // Try HEAD first (lighter), fall back to GET
    let init = RequestInit::new();
    init.set_method("HEAD");
    init.set_mode(RequestMode::Cors);

    let response_promise: js_sys::Promise = window.fetch_with_str_and_init(url, &init).into();
    let response_value = match JsFuture::from(response_promise).await
    {
        Ok(v) => v,
        Err(_) => return false,
    };

    if let Ok(response) = response_value.dyn_into::<Response>()
    {
        response.ok()
    }
    else
    {
        false
    }
}

/// Fetch both index files for a FASTA URL
/// Returns (fai_data, gzi_data) - FAI is required, GZI is optional
/// Succeeds if at least FAI can be fetched
#[cfg(target_arch = "wasm32")]
pub async fn fetch_index_files(fasta_url: &str) -> Result<(Vec<u8>, Option<Vec<u8>>)>
{
    let fai_url = format!("{}.fai", fasta_url);
    let gzi_url = format!("{}.gzi", fasta_url);

    // Fetch FAI (required)
    let fai_data = load_file_async(&fai_url).await?;

    // Fetch GZI (only for gzipped files, optional but recommended)
    let gzi_data = if fasta_url.ends_with(".gz")
    {
        match load_file_async(&gzi_url).await
        {
            Ok(data) => Some(data),
            Err(_) => None, // GZI not available, but we can still work with just FAI for uncompressed
        }
    }
    else
    {
        None
    };

    Ok((fai_data, gzi_data))
}
