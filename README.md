# ARPES.jl

[![Build Status](https://github.com/arafune/ARPES.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/arafune/ARPES.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/arafune/ARPES.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/arafune/ARPES.jl)

ARPES.jl is a Julia package for working with **Angle-Resolved Photoemission Spectroscopy (ARPES)** data.
It provides core data structures, dimensional axes utilities, I/O helpers, and **k-space conversion** utilities.
It aims to enable consistent workflows via `ARPESData`, an extension of [DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl) based data container.

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

# Load ARPES data (file formats depend on the IO backends implemented in `src/io/`)
data = load("path/to/data", loc="LOCATION")

# k-space conversion (see `k_conversion`)
kdata = k_conversion(data)
```

---

## Public API (exported)

From `src/ARPES.jl`, the package exports:

- Axes / dimensions: `kx`, `ky`, `kz`, `phi`, `psi`, `eV`, `detector_ch`, `ch2`, `delay`, `spin`
  - `phi` : Emission angle of photoelectron along the analyzer slit direction.
  - `psi` : Emission angle of photoelectron perpendicular to the analyzer slit direction.
  - `kx`, `ky`, `kz` : Momentum components of the photoelectron, defined in the sample's reciprocal space.
  - `eV` : Energy of the photoelectron in electron volts.

- Units / signals: `CPS`, `Counts`
- Core data type: `ARPESData`
- I/O: `load`
- k-space conversion: `k_conversion`

---

## Repository layout

- `src/` — core implementation
  - `ARPES.jl` — main module and exports
  - `types.jl` — core type definitions
  - `ARPESData.jl` — ARPES data container
  - `dims.jl` — dimension/axis helpers
  - `io/` — loaders (`load`)
  - `k_conv/` — k-space conversion (`k_conversion`)
- `docs/` — documentation (Documenter.jl)
- `test/` — tests
- `testdata` — git submodule (test fixtures)

---

## Development with Dev Containers

This repository includes a Dev Container configuration for easy development setup. You can use it with:

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

## Development with Dev Containers

This repository includes a Dev Container configuration for easy development setup. You can use it with:

- **VS Code + Docker**: Open in VS Code and click "Reopen in Container" when prompted
- **GitHub Codespaces**: Click the "Code" button on GitHub and create a Codespace

The Dev Container automatically sets up:

- Julia 1.12 environment
- All package dependencies
- VS Code extensions for Julia development
- Test and documentation environments

For detailed instructions, see [.devcontainer/README.md](.devcontainer/README.md).

[![codecov](https://codecov.io/gh/arafune/ARPES.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/arafune/ARPES.jl)
