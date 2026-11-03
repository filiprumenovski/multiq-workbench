//! Colours, fonts and egui visuals: dark ink by default, warm paper in light mode, one accent hue, muted species colours.
use egui::{self, Color32, FontData, FontDefinitions, FontFamily, FontId, Stroke, TextStyle};
use std::sync::Arc;

#[derive(Clone, Copy)]
pub struct Palette {
    pub dark: bool,
    pub bg: Color32,       // canvas
    pub panel: Color32,    // side panels
    pub card: Color32,     // raised rows
    pub line: Color32,     // hairlines / dividers
    pub text: Color32,
    pub dim: Color32,
    pub faint: Color32,
    pub contour: Color32,
    pub index: Color32,
    pub elec_fill: Color32,
    pub elec_line: Color32,
    pub ghost: Color32,
    pub accent: Color32,
    pub ok: Color32,       // settled / consistent
    pub bad: Color32,      // drifting / inconsistent
}

pub fn palette(dark: bool) -> Palette {
    if dark {
        Palette {
            dark, bg: rgb(0x0B0D10), panel: rgb(0x0F1216), card: rgb(0x151A20), line: rgb(0x242B34),
            text: rgb(0xD8DEE6), dim: rgb(0x8B96A5), faint: rgb(0x4A5460),
            contour: Color32::from_rgba_unmultiplied(0x6F, 0x86, 0x9C, 0x62), index: Color32::from_rgba_unmultiplied(0xA6, 0xB8, 0xCB, 0xC8),
            elec_fill: rgb(0x1B2128), elec_line: rgb(0x7D8A99), ghost: Color32::from_rgba_unmultiplied(0x7D, 0x8A, 0x99, 0x38),
            accent: rgb(0xF2B04A), ok: rgb(0x6FBF8E), bad: rgb(0xD97B7B),
        }
    } else {
        Palette {
            dark, bg: rgb(0xF4F1EA), panel: rgb(0xEDE9E0), card: rgb(0xE5E0D5), line: rgb(0xD0CABD),
            text: rgb(0x23272D), dim: rgb(0x68707B), faint: rgb(0xA6A79F),
            contour: Color32::from_rgba_unmultiplied(0x6A, 0x5A, 0x44, 0x6A), index: Color32::from_rgba_unmultiplied(0x4A, 0x3C, 0x2A, 0xE0),
            elec_fill: rgb(0xCFC9BB), elec_line: rgb(0x5A5348), ghost: Color32::from_rgba_unmultiplied(0x5A, 0x53, 0x48, 0x34),
            accent: rgb(0xC9731A), ok: rgb(0x2F8F5F), bad: rgb(0xB0484A),
        }
    }
}

const fn rgb(c: u32) -> Color32 { Color32::from_rgb((c >> 16) as u8, (c >> 8) as u8, c as u8) }

/// Muted, distinguishable species colours (cool family; the accent amber is reserved for excitation).
pub fn species_color(dark: bool, i: usize) -> Color32 {
    const D: [u32; 6] = [0x56B6A6, 0x6E9BDA, 0xA68BDA, 0xD27FB0, 0x8FC27A, 0x58B8D6];
    const L: [u32; 6] = [0x2F8F80, 0x3F6FB5, 0x7A5FB8, 0xB0568C, 0x5E9A45, 0x2D8CAD];
    let c = if dark { D[i % 6] } else { L[i % 6] };
    rgb(c)
}

pub fn with_alpha(c: Color32, a: u8) -> Color32 { Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), a) }

pub const MONO: &str = "mq-mono";
pub const SANS: &str = "mq-sans";

pub fn mono(size: f32) -> FontId { FontId::new(size, FontFamily::Name(MONO.into())) }
pub fn sans(size: f32) -> FontId { FontId::new(size, FontFamily::Name(SANS.into())) }

fn load(fonts: &mut FontDefinitions, name: &str, path: &str, fallback: FontFamily) -> bool {
    if let Ok(bytes) = std::fs::read(path) {
        fonts.font_data.insert(name.to_owned(), Arc::new(FontData::from_owned(bytes)));
        // keep egui's built-in glyph fallbacks behind the system font
        let mut list = vec![name.to_owned()];
        list.extend(fonts.families.get(&fallback).cloned().unwrap_or_default());
        fonts.families.insert(FontFamily::Name(name.into()), list);
        true
    } else { false }
}

pub fn setup_fonts(ctx: &egui::Context) {
    let mut fonts = FontDefinitions::default();
    // one mono font for data (SF Mono -> Menlo-like fallback to egui's Hack), one clean sans for labels
    let mono_ok = load(&mut fonts, MONO, "/System/Library/Fonts/SFNSMono.ttf", FontFamily::Monospace);
    let sans_ok = load(&mut fonts, SANS, "/System/Library/Fonts/SFNS.ttf", FontFamily::Proportional);
    if !mono_ok { let fb = fonts.families.get(&FontFamily::Monospace).cloned().unwrap_or_default(); fonts.families.insert(FontFamily::Name(MONO.into()), fb); }
    if !sans_ok { let fb = fonts.families.get(&FontFamily::Proportional).cloned().unwrap_or_default(); fonts.families.insert(FontFamily::Name(SANS.into()), fb); }
    ctx.set_fonts(fonts);
}

pub fn apply(ctx: &egui::Context, p: &Palette) {
    let mut v = if p.dark { egui::Visuals::dark() } else { egui::Visuals::light() };
    v.panel_fill = p.panel;
    v.window_fill = p.panel;
    v.extreme_bg_color = p.bg;
    v.faint_bg_color = p.card;
    v.override_text_color = Some(p.text);
    v.window_stroke = Stroke::new(1.0, p.line);
    v.selection.bg_fill = with_alpha(p.accent, 90);
    v.selection.stroke = Stroke::new(1.0, p.accent);
    v.hyperlink_color = p.accent;
    let w = &mut v.widgets;
    let rail = if p.dark { rgb(0x252D37) } else { rgb(0xD8D2C4) };
    for (s, bg) in [(&mut w.inactive, rail), (&mut w.hovered, p.line), (&mut w.active, p.line), (&mut w.open, rail)] {
        s.bg_fill = bg;
        s.weak_bg_fill = bg;
        s.bg_stroke = Stroke::new(1.0, p.line);
        s.fg_stroke = Stroke::new(1.0, p.text);
        s.corner_radius = egui::CornerRadius::same(3);
    }
    w.noninteractive.bg_stroke = Stroke::new(1.0, p.line);
    w.noninteractive.fg_stroke = Stroke::new(1.0, p.dim);
    w.hovered.bg_stroke = Stroke::new(1.0, p.dim);
    w.active.bg_stroke = Stroke::new(1.0, p.accent);
    ctx.set_visuals(v);
    ctx.global_style_mut(|st| {
        st.text_styles.insert(TextStyle::Body, sans(12.5));
        st.text_styles.insert(TextStyle::Button, sans(12.5));
        st.text_styles.insert(TextStyle::Small, sans(10.5));
        st.text_styles.insert(TextStyle::Monospace, mono(11.5));
        st.text_styles.insert(TextStyle::Heading, sans(15.0));
        st.spacing.item_spacing = egui::vec2(8.0, 5.0);
        st.spacing.slider_width = 150.0;
        st.spacing.interact_size.y = 20.0;
    });
}
