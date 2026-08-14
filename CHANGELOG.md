# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.3] — Unreleased

### Added

- `tarpes_evolution` (src/pumpprobe.jl): helper to extract a 2D snapshot and the temporal evolution from 3D time-resolved ARPES datasets. Supports selection by delay time or delay index, scalar or (center, width) averaging along the non-dispersion axis, and a `full_temporal` flag to include the full time range. Unit tests added in `test/pumpprobe.jl`.

### Changed

- Update dimensions, preprocessing, and trapezoid calculations (updates in src/dims.jl, src/k_conv/preprocess.jl, and src/trapezoid.jl) (#71).
- CI: bump actions and workflows (setup-julia → v3, codecov action → v7, actions/checkout → v7) (#68, #69, #70).

### Fixed

- Refactored src/io/formats/itx.jl and src/k_conv/preprocess.jl: reorganized ITX parsing into clearer helper functions and improved robustness of scale/waves parsing and DimArray construction; added helpers in k_conv/preprocess.jl (prepare_for_broadcast, _ek_range, _kx_range, _ky_range) and stricter ARPESData validation (_check_arpesdata) to make k-conversion preprocessing more deterministic. Corresponding tests added/updated in test/io/formats/itx.jl to cover WAVES/SetScale/comment parsing and wave data assembly (commit 36aa831).

## [0.2.2] — 2026-04-16

### Added

- `dim_extend` keyword argument for `shift`
  - Supports scalar, vector, and one-dimensional `AbstractDimArray` shift inputs.
  - Extends equally spaced shift dimensions so shifted data is preserved on a larger output axis.
- `trapezoid` for 2D `ARPESData` with `:phi` and `:eV` dimensions ([#66]).
  - Supports both rectangle → trapezoid and trapezoid → rectangle coordinate correction.
  - Uses linear interpolation on angular slices and preserves the `:eV` axis.
  - Added unit tests for coordinate mapping, both conversion directions, and invalid geometry handling.

### Changed

- Consolidated and expanded the `shift` docstrings to document scalar shifts, per-dimension shifts, irregular-grid behavior, and `dim_extend`.

## [0.2.1] — 2026-04-13

### Added

- `stitch_along`
  - Stitch two `AbstractDimArray`s together along a specified dimension.
  - Handles overlapping regions: `seam_ratio` controls the split point within the overlap (0–1; `nothing` concatenates directly without splitting).
  - `gain_a` applies a multiplicative scale to the first array before concatenating.
  - Result is always sorted along the stitched dimension.

### Fixed

- Fixed insufficient parse option in SetScale in `load_itx` ([#63]).

## [0.2.0] — 2026-04-08

### Breaking change (Removed)

- `cat_arpes` removed; the standard `Base.cat` is sufficient ([#52]).

### Added

- **3D `k_conversion`** for `ARPESData{T,3}` ([#52]):
  - Dispatches to `kxky_conversion` when a `:psi` dimension is present (kx-ky map).
  - Dispatches to `kpkz_conversion` when a `:hv`/`:hν` dimension is present (kp-kz map, stub).
  - Otherwise iterates over the third scan dimension, applies the 2D `k_conversion` on each φ-eV slice, and concatenates the results.
- **4D `k_conversion`** for `ARPESData{T,4}` ([#52]):
  - Dispatches to `kxkykz_conversion` (stub) for φ+ψ+hν data.
  - Falls back to iterating over the extra scan dimension and delegating to 3D dispatch.
- **N≥5D `k_conversion`** ([#52]):
  - Recursively peels the last scan dimension (any dimension that is not eV, phi, psi, hv, or hν) and delegates to the (N-1)D dispatch.
- `kxky_conversion` helper for 3D φ×ψ×eV → kx×ky×eV conversion ([#52]).
- Stubs for `kpkz_conversion` and `kxkykz_conversion` ([#52]).
- `_slice_and_convert` internal helper to eliminate repeated slice-stack patterns across dispatch levels ([#52]).
- 4D quadrilinear interpolation overloads for `_interpolate` (uniform and non-uniform grids) ([#52]).
- `ky_range` and `kz_range` keyword arguments added to the 2D `k_conversion` signature for a consistent call interface from higher-dimensional dispatches ([#52]).
- `AbstractDimArray` and `DimensionalData.Dimension` overloads for `negate_dim` and `add_dim` in `dims.jl` ([#52]).
- `testdata` added as a git submodule (dev branch only) ([#50]).

### Changed

- Renamed `reshape_for_nd` → `prepare_for_broadcast` ([#52]).
- Renamed loop variable `a_dimarray` → `arpes_cut` in kx-ky dispatch for clarity ([#52]).
- Renamed `_interpolate` parameters from 0-based (`p0`/`c0`) to 1-based (`p1`/`c1`) to align with Julia's indexing convention ([#52]).
- Renamed local variable `dims` → `dim_labels` in `shift_dim` varargs to avoid shadowing ([#52]).
- Renamed `otherdim` → `target_dim` throughout N=3 and N=4 fallback cases in `k_conversion` ([#52]).
- `k_conversion` docstring rewritten with a full dispatch table covering N=2 through N≥5 cases ([#52]).
- `derivative.jl`: added missing `return` in `curvature/2` (previously returned `nothing`) ([#52]).
- `filter.jl`: hoisted `weights` allocation out of the inner loop in `sg_1d`; use `@view` for slices to avoid copies ([#52]).
- CI: added `deploy-main.yml` workflow to maintain a clean `main` branch without testdata ([#50]).
- Tests split into two sections for cleaner CI reporting ([#54]).

### Fixed

- Fixed loading of n-dimensional (n > 2) ITX data ([#45]).
- Fixed `ky_range` being always overwritten in 2D `k_conversion` (added `isnothing` guard) ([#52]).
- Fixed operator-precedence bug in `hasdim(:hv) || hasdim(:hν)` check ([#52]).
- Fixed `UndefVarError` in `add_dim` caused by misnamed parameter (`dim_label` → `dim`) ([#52]).
- Fixed `UndefVarError` references to undefined `otherdim` in N=3 and N=4 fallback paths ([#52]).
- Fixed dimension-order reversal in the N=4 fallback: changed `others[1]` to `others[end]` so recursion preserves input order (e.g. φ, eV, A, B → kx, eV, A, B) ([#52]).
- Fixed `_shift_index` docstring header (was incorrectly labeled `_process_index`) ([#52]).
- Fixed typos

---

## [0.1.1] — 2026-04-02

### Added

- `rebin` function with documentation ([#38]).
- Positional overloads and defaults for filter functions ([#37-area]).
- `convert_dim` for `DimensionalData.Dimension` including overloads for default-label and `String`-label usage.
- Improved handling of the angle χ in k-conversion (`_set_supplemental_angles`).
- GitHub Pages documentation deployment.
- Expanded documentation guide and API reference.

### Changed

- Simplified k-conversion dimension checks.
- Simplified `convert_dim` API.
- Refactored `AnalyzerConfiguration` type hierarchy (new abstract type below `AnalyzerConfiguration`).
- Various code style and refactoring improvements.

### Fixed

- CompatHelper workflow corrections ([#37]).
- Removed `LinearAlgebra` and `Statistics` from `[deps]` to fix CompatHelper.

[0.1.1]: https://github.com/arafune/ARPES.jl/releases/tag/v0.1.1
[#45]: https://github.com/arafune/ARPES.jl/pull/45
[#50]: https://github.com/arafune/ARPES.jl/pull/50
[#52]: https://github.com/arafune/ARPES.jl/pull/52
[#54]: https://github.com/arafune/ARPES.jl/pull/54
[#66]: https://github.com/arafune/ARPES.jl/pull/66
