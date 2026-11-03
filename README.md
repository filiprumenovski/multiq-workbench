# MultiQ Workbench

A macOS desktop simulator for the [MultiQ ion trap](https://doi.org/10.1126/sciadv.aec7048). Explore shows computed ion trajectories through one channel. Tandem tracks a shared ion population through a splitter trap, nine channels, and a combiner trap. An ion rejected by one channel can return to the trap and enter another.

## Run

Install the Xcode Command Line Tools, stable Rust, Node.js with npm, and the Tauri 2 CLI. From the repository root:

```sh
cargo install tauri-cli --version '^2' --locked  # once
npm ci
npm run dev
```

Choose **Tandem → Run tandem** to start with four channels driven for the Fig. 5 ubiquitin experiment. Set a channel’s frequency to zero to block it. Save a run to export the settings, charge-state counts and port-to-port retries. Use `npm run bundle` to build the macOS app.

## What has been checked

- **Original code:** the [released velocity update](https://github.com/akrutchins/Ion-Ejection-from-MultiQ/blob/279240a314978ea508aa43e60c8bd212b11321d8/TRAJECTORY.pde#L90-L97) advances velocity by half the acceleration per step. At m/z 499, both original Java and legacy Rust reproduced 742 transmissions per 1,000 launches. With the full velocity update, neither transmitted an ion anywhere in the tested m/z 475–523 scan at those drive settings.
- **Numerical sensitivity:** on the same 1,000 m/z 499 launches, the corrected model returned 175, 196 and 207 ions on fresh 0.10, 0.05 and 0.0333 mm force grids, respectively. Most ions hit an electrode at all three spacings; the return fraction is not converged.

The bundled data are digitized published figures. The model has not been validated against raw instrument scans or transients. The trap cloud in the instrument render is illustrative. Tandem uses independent ions and separately solved trap and channel fields; it omits space charge. The next experimental check needs raw scans and transients with their settings and timing.

## Tests and sources

```sh
npm run build && npm test
cargo test --release --locked --manifest-path src-tauri/Cargo.toml -p multiq-desktop -p trapsim
```

The original [486Q trap](https://github.com/akrutchins/MultiQ) and [ejection channel](https://github.com/akrutchins/Ion-Ejection-from-MultiQ) sketches by Andrew Krutchinsky are in [reference/](reference/) with their BSD notices. Opening those sketches uses Processing 4; the Rust simulator runs without it.
