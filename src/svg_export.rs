use crate::renderer::TrackPainter;
use egui::{Align2, Color32, FontFamily, FontId, Pos2, Rect, Stroke, Vec2};
use std::cell::RefCell;

pub struct SvgPainter
{
    origin: Pos2,
    size: Vec2,
    elements: RefCell<Vec<String>>,
}

impl SvgPainter
{
    pub fn new(rect: Rect) -> Self
    {
        Self {
            origin: rect.min,
            size: rect.size(),
            elements: RefCell::new(Vec::new()),
        }
    }

    pub fn finish(self) -> String
    {
        let mut svg = String::new();
        svg.push_str(&format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">"#,
            w = self.size.x.max(1.0),
            h = self.size.y.max(1.0)
        ));
        for element in self.elements.into_inner()
        {
            svg.push('\n');
            svg.push_str(&element);
        }
        svg.push_str("\n</svg>");
        svg
    }

    fn push(&self, element: String)
    {
        self.elements.borrow_mut().push(element);
    }

    fn to_local(&self, pos: Pos2) -> Pos2
    {
        Pos2::new(pos.x - self.origin.x, pos.y - self.origin.y)
    }

    fn escape_text(text: &str) -> String
    {
        text.replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
            .replace('"', "&quot;")
            .replace('\'', "&apos;")
    }

    fn color_to_rgb(color: Color32) -> String
    {
        format!("#{:02x}{:02x}{:02x}", color.r(), color.g(), color.b())
    }

    fn color_opacity(color: Color32) -> f32
    {
        (color.a() as f32 / 255.0).clamp(0.0, 1.0)
    }

    fn align_to_svg(align: Align2) -> (&'static str, &'static str)
    {
        let anchor = match align.x()
        {
            egui::Align::Min => "start",
            egui::Align::Center => "middle",
            egui::Align::Max => "end",
        };
        let baseline = match align.y()
        {
            egui::Align::Min => "hanging",
            egui::Align::Center => "middle",
            egui::Align::Max => "text-after-edge",
        };
        (anchor, baseline)
    }
}

impl TrackPainter for SvgPainter
{
    fn rect_filled(&self, rect: Rect, rounding: f32, color: Color32)
    {
        let min = self.to_local(rect.min);
        let size = rect.size();
        let opacity = Self::color_opacity(color);
        let mut element = format!(
            r#"<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" ry="{r}" fill="{fill}""#,
            x = min.x,
            y = min.y,
            w = size.x.max(0.0),
            h = size.y.max(0.0),
            r = rounding.max(0.0),
            fill = Self::color_to_rgb(color)
        );
        if opacity < 1.0
        {
            element.push_str(&format!(r#" fill-opacity="{opacity}""#, opacity = opacity));
        }
        element.push_str(" />");
        self.push(element);
    }

    fn rect_stroke(&self, rect: Rect, rounding: f32, stroke: Stroke)
    {
        let min = self.to_local(rect.min);
        let size = rect.size();
        let opacity = Self::color_opacity(stroke.color);
        let mut element = format!(
            r#"<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" ry="{r}" fill="none" stroke="{stroke}" stroke-width="{sw}""#,
            x = min.x,
            y = min.y,
            w = size.x.max(0.0),
            h = size.y.max(0.0),
            r = rounding.max(0.0),
            stroke = Self::color_to_rgb(stroke.color),
            sw = stroke.width
        );
        if opacity < 1.0
        {
            element.push_str(&format!(r#" stroke-opacity="{opacity}""#, opacity = opacity));
        }
        element.push_str(" />");
        self.push(element);
    }

    fn line_segment(&self, points: [Pos2; 2], stroke: Stroke)
    {
        let p0 = self.to_local(points[0]);
        let p1 = self.to_local(points[1]);
        let opacity = Self::color_opacity(stroke.color);
        let mut element = format!(
            r#"<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{sw}""#,
            x1 = p0.x,
            y1 = p0.y,
            x2 = p1.x,
            y2 = p1.y,
            stroke = Self::color_to_rgb(stroke.color),
            sw = stroke.width
        );
        if opacity < 1.0
        {
            element.push_str(&format!(r#" stroke-opacity="{opacity}""#, opacity = opacity));
        }
        element.push_str(" />");
        self.push(element);
    }

    fn text(&self, pos: Pos2, align: Align2, text: &str, font: FontId, color: Color32)
    {
        let p = self.to_local(pos);
        let (anchor, baseline) = Self::align_to_svg(align);
        let opacity = Self::color_opacity(color);
        let family = match font.family
        {
            FontFamily::Proportional => "sans-serif",
            FontFamily::Monospace => "monospace",
            FontFamily::Name(ref name) => name.as_ref(),
        };
        let mut element = format!(
            r#"<text x="{x}" y="{y}" font-size="{size}" font-family="{family}" text-anchor="{anchor}" dominant-baseline="{baseline}" fill="{fill}""#,
            x = p.x,
            y = p.y,
            size = font.size,
            family = family,
            anchor = anchor,
            baseline = baseline,
            fill = Self::color_to_rgb(color)
        );
        if opacity < 1.0
        {
            element.push_str(&format!(r#" fill-opacity="{opacity}""#, opacity = opacity));
        }
        element.push_str(">");
        element.push_str(&Self::escape_text(text));
        element.push_str("</text>");
        self.push(element);
    }

    fn circle_filled(&self, center: Pos2, radius: f32, color: Color32)
    {
        let c = self.to_local(center);
        let opacity = Self::color_opacity(color);
        let mut element = format!(
            r#"<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}""#,
            x = c.x,
            y = c.y,
            r = radius.max(0.0),
            fill = Self::color_to_rgb(color)
        );
        if opacity < 1.0
        {
            element.push_str(&format!(r#" fill-opacity="{opacity}""#, opacity = opacity));
        }
        element.push_str(" />");
        self.push(element);
    }
}
