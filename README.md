# ARPES.jl

[![Build Status](https://github.com/arafune/ARPES.jl/actions/workflows/CI.yml/badge.svg?branch=dev)](https://github.com/arafune/ARPES.jl/actions/workflows/CI.yml?query=branch%3Adev)
[![codecov](https://codecov.io/gh/arafune/ARPES.jl/branch/dev/graph/badge.svg)](https://codecov.io/gh/arafune/ARPES.jl)

ARPES.jl is a Julia package for working with
**Angle Resolved Photoemission Spectroscopy (ARPES)** data.
It provides a dimension-aware `ARPESData` container built
on top of [DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl),
together with utilities for loading, transforming axes, filtering,
differentiating, and converting ARPES datasets into momentum space.

Key capabilities in the current `src/` tree:

- canonical ARPES dimensions such as `phi`, `psi`, `eV`, `kx`, `ky`, and `kz`
- `ARPESData`, a thin wrapper around `AbstractDimArray` that preserves
  dimensional indexing and metadata
- a layered I/O pipeline for parsing raw files and standardizing them into the
  package's canonical conventions
- dimension and coordinate helpers for renaming, shifting, and transforming axes
- signal-processing helpers for smoothing, derivatives, and curvature-style
  feature enhancement
- `k_conversion` for converting angle-energy data into momentum-space coordinates

> Julia compatibility: Julia 1.12 (see `Project.toml`)

---

## Installation

```julia
using Pkg
Pkg.add("ARPES")
```

or for development:

```julia
using Pkg
Pkg.develop(url="https://github.com/arafune/ARPES.jl")
```

---

## Quick start

```julia
using ARPES

# Load data with the currently implemented SPD loader
data = load("path/to/data.itx"; loc = "SPD")

# Dimension-aware utilities
shifted = shift_dim(data, :phi, 0.2)
rescaled = convert_dim(data, :phi, x -> 0.5x)
dIdE = derivative(data, :eV)
curv = curvature(data, :eV)

# k-space conversion
kdata = k_conversion(data)
```

---

## Current I/O support

The current I/O implementation is centered on the `SPDLoader`.

- Supported file types: `.itx`, `.sp2`
- Loader selection:
  - pass `loc = "SPD"` explicitly, or
  - omit `loc` and let the registry detect the loader from the file extension
    and ITX header
- Loading pipeline:
  - `src/io/formats/` parses the raw file into a mostly faithful `DimArray`
  - `src/io/location/` maps location-specific axes and metadata onto canonical
    ARPES names
  - the standardized result is wrapped as `ARPESData`

At the moment, SPD is the concrete location-specific loader registered in `src/io/registry.jl`.

---

## Public API

The package currently exports the following user-facing names:

- Dimensions: `kx`, `ky`, `kz`, `phi`, `psi`, `eV`, `detector_ch`, `ch2`,
  `delay`, `spin`
- Intensity units: `CPS`, `Counts`
- Core type: `ARPESData`
- Metadata helper: `merge_consensus`
- Dimension helpers: `add_dim`, `convert_dim`, `shift_dim`, `negate_dim`, `rename_dim`
- Data transforms: `shift`, `rebuild_with_slice`
- Filters: `kalman_smooth_dim_llpf`, `boxcar_filter_dim`, `sg_filter_dim`,
  `binomial_filter_dim`, `gaussian_filter_dim`
- Derivative / feature extraction: `derivative`, `curvature`, `minimum_gradient`
- I/O: `load`
- Momentum conversion: `k_conversion`

Some important conventions used across the package:

- `phi` is the emission angle along the analyzer slit direction
- `psi` is the emission angle perpendicular to the slit direction
- `eV` is the energy axis, with the actual energy meaning carried by
  `dims(data, :eV)` metadata via `:energy_definition`
- `convert_dim` preserves the dimension label and metadata while transforming
  coordinate values
- dataset-level metadata is part of the package contract and is used by loaders
  and `k_conversion`

---

## Repository layout

- `src/` — core implementation
  - `ARPES.jl` — main module and exports
  - `types.jl` — core type definitions
  - `ARPESData.jl` — ARPES data container
  - `dims.jl` — dimension/axis helpers
  - `transform.jl` — shift and slice-based data transforms
  - `filter.jl` — smoothing utilities along selected dimensions
  - `derivative.jl` — numerical derivatives and 2D feature-enhancement helpers
  - `io/` — loaders (`load`)
    - `formats/` — raw file parsing
    - `location/` — location-specific standardization into canonical ARPES conventions
  - `k_conv/` — k-space conversion (`k_conversion`)
- `docs/` — documentation (Documenter.jl)
- `test/` — tests

---

## Development

Instantiate and build the package environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
```

Run the full test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run a single test file directly from the repository root:

```bash
julia --project=. -e 'include("test/io/io.jl")'
```

Build the documentation:

```bash
julia --project=docs docs/make.jl
```

Format the source tree with the repository's JuliaFormatter configuration:

```bash
julia --project=. -e 'using JuliaFormatter; format("src")'
```

### Development with Dev Containers

This repository includes a Dev Container configuration for easy development
setup. You can use it with:

- **VS Code + Docker**: Open in VS Code and click "Reopen in Container" when prompted
- **GitHub Codespaces**: Click the "Code" button on GitHub and create a Codespace

The Dev Container automatically sets up:

- Julia 1.12 environment
- All package dependencies
- VS Code extensions for Julia development
- Test and documentation environments

For detailed instructions, see [.devcontainer/README.md](.devcontainer/README.md).

---

## License

See [LICENSE](LICENSE).
