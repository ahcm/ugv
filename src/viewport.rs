#[derive(Debug, Clone)]
pub struct Viewport
{
    pub start: usize,
    pub end: usize,
    max_length: usize,
}

impl Viewport
{
    pub fn new(start: usize, max_length: usize) -> Self
    {
        let initial_width = max_length.min(1_000_000);
        Self {
            start,
            end: start + initial_width,
            max_length,
        }
    }

    pub fn width(&self) -> usize
    {
        self.end.saturating_sub(self.start)
    }

    pub fn pan(&mut self, delta: i64)
    {
        let new_start = (self.start as i64 + delta).max(0) as usize;
        let width = self.width();

        if new_start + width <= self.max_length
        {
            self.start = new_start;
            self.end = new_start + width;
        }
        else
        {
            self.end = self.max_length;
            self.start = self.max_length.saturating_sub(width);
        }
    }

    pub fn zoom(&mut self, factor: f64, focus: f32)
    {
        let width = self.width() as f64;
        let new_width = (width * factor).max(100.0).min(self.max_length as f64) as usize;

        let focus_point = self.start + (width * focus as f64) as usize;
        let new_start = focus_point.saturating_sub((new_width as f64 * focus as f64) as usize);

        if new_start + new_width <= self.max_length
        {
            self.start = new_start;
            self.end = new_start + new_width;
        }
        else
        {
            self.end = self.max_length;
            self.start = self.max_length.saturating_sub(new_width);
        }
    }

    pub fn set_region(&mut self, start: usize, end: usize)
    {
        self.start = start.min(self.max_length);
        self.end = end.min(self.max_length);
    }

    pub fn zoom_to_feature(&mut self, feature_start: usize, feature_end: usize, padding: f32)
    {
        let feature_width = feature_end.saturating_sub(feature_start);
        let padded_width = (feature_width as f32 * (1.0 + padding * 2.0)) as usize;

        let center = (feature_start + feature_end) / 2;
        let new_start = center.saturating_sub(padded_width / 2);
        let new_end = (new_start + padded_width).min(self.max_length);

        self.start = new_start;
        self.end = new_end;
    }

    pub fn position_to_screen(&self, position: usize, screen_width: f32) -> f32
    {
        let relative_pos = position.saturating_sub(self.start) as f32;
        let viewport_width = self.width() as f32;
        (relative_pos / viewport_width) * screen_width
    }

    pub fn screen_to_position(&self, screen_x: f32, screen_width: f32) -> usize
    {
        let ratio = screen_x / screen_width;
        self.start + (self.width() as f32 * ratio) as usize
    }
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn test_viewport_pan()
    {
        let mut vp = Viewport::new(0, 10000);
        vp.end = 1000;

        vp.pan(500);
        assert_eq!(vp.start, 500);
        assert_eq!(vp.end, 1500);

        vp.pan(-200);
        assert_eq!(vp.start, 300);
        assert_eq!(vp.end, 1300);
    }

    #[test]
    fn test_viewport_zoom()
    {
        let mut vp = Viewport::new(0, 10000);
        vp.end = 1000;

        vp.zoom(0.5, 0.5);
        assert!(vp.width() < 1000);

        vp.zoom(2.0, 0.5);
        assert!(vp.width() > 500);
    }
}
