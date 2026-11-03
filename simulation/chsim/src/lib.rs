//! chsim as a library: field solves (`channel::Fields`, `channel::Geom`, `s14_geom`) and the Metal engine, so that
//! other crates (mqviz) can reuse the validated channel geometry and field code. The `chsim` binary keeps its own
//! module tree (src/main.rs) and is unaffected.
pub mod channel;
pub mod check;
pub mod tiled;
// spectrum forward model (SPECTRUM.md): CPU-only, reused by mqviz's reality panel (predicted single-channel spectra)
pub mod scan;
pub mod spectrum;
pub mod spectrum_detector;
pub mod spectrum_like;
pub mod spectrum_trans;
pub mod spectrum_val;
