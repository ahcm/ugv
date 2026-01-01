/// Generate the application icon (DNA helix design)
pub fn create_app_icon() -> egui::IconData
{
    let size = 64;
    let mut rgba = vec![0u8; size * size * 4];

    // DNA helix colors
    let bg_color = [245, 247, 250, 255]; // Light gray background
    let strand1_color = [70, 130, 180, 255]; // Steel blue
    let strand2_color = [220, 100, 100, 255]; // Coral red
    let connector_color = [150, 150, 150, 255]; // Gray

    // Fill background
    for y in 0..size
    {
        for x in 0..size
        {
            let idx = (y * size + x) * 4;
            rgba[idx..idx + 4].copy_from_slice(&bg_color);
        }
    }

    // Draw DNA double helix
    let center_x = size as f32 / 2.0;
    let amplitude = 18.0;
    let frequency = 0.15;

    for y in 0..size
    {
        let t = y as f32 * frequency;

        // First strand (sine wave)
        let x1 = center_x + amplitude * t.sin();
        // Second strand (cosine wave - 180° out of phase)
        let x2 = center_x - amplitude * t.sin();

        // Draw first strand
        draw_circle(&mut rgba, size, x1 as i32, y as i32, 3, &strand1_color);

        // Draw second strand
        draw_circle(&mut rgba, size, x2 as i32, y as i32, 3, &strand2_color);

        // Draw connectors at crossover points (every ~20 pixels)
        if y % 10 == 0 && y > 0 && y < size - 1
        {
            draw_line(
                &mut rgba,
                size,
                x1 as i32,
                y as i32,
                x2 as i32,
                y as i32,
                &connector_color,
            );
        }
    }

    egui::IconData {
        rgba,
        width: size as u32,
        height: size as u32,
    }
}

fn draw_circle(rgba: &mut [u8], size: usize, cx: i32, cy: i32, radius: i32, color: &[u8; 4])
{
    for dy in -radius..=radius
    {
        for dx in -radius..=radius
        {
            if dx * dx + dy * dy <= radius * radius
            {
                let x = cx + dx;
                let y = cy + dy;

                if x >= 0 && x < size as i32 && y >= 0 && y < size as i32
                {
                    let idx = (y as usize * size + x as usize) * 4;
                    rgba[idx..idx + 4].copy_from_slice(color);
                }
            }
        }
    }
}

fn draw_line(
    rgba: &mut [u8],
    size: usize,
    x0: i32,
    y0: i32,
    x1: i32,
    y1: i32,
    color: &[u8; 4],
)
{
    let dx = (x1 - x0).abs();
    let dy = (y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx - dy;

    let mut x = x0;
    let mut y = y0;

    loop
    {
        if x >= 0 && x < size as i32 && y >= 0 && y < size as i32
        {
            let idx = (y as usize * size + x as usize) * 4;
            rgba[idx..idx + 4].copy_from_slice(color);
        }

        if x == x1 && y == y1
        {
            break;
        }

        let e2 = 2 * err;
        if e2 > -dy
        {
            err -= dy;
            x += sx;
        }
        if e2 < dx
        {
            err += dx;
            y += sy;
        }
    }
}
