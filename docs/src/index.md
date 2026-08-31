```@meta
CurrentModule = ARPES
```

# ARPES.jl

`ARPES.jl` is a Julia package for working with angle-resolved photoemission spectroscopy
(ARPES) data on
[DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl).

The package is built around three ideas:

- keep spectra as dimension-aware arrays rather than plain matrices,
- preserve experimental metadata as part of the data contract, and
- provide a direct path from raw angle-space data to momentum-space data.

## What the package provides

- `ARPESData` as the main container for spectra and metadata,
- canonical dimensions such as `eV`, `phi`, `psi`, `kx`, `ky`, and `kz`,
- file loading through `load`,
- dimension-aware axis editing, transforms, smoothing, derivatives, and curvature
  enhancement, and
- `k_conversion` for angle-to-momentum conversion.

## What the package does _not_ provide

- Plotting utilities, which are best left to dedicated packages such as
  [Makie.jl](https://makie.juliaplots.org/stable/). See [ARPESMakie.jl](https://arafune.github.io/ARPESMakie.jl).

## Quick start

```julia
using ARPES

# Load a dataset. The current concrete location-specific loader is SPD.
data = load("path/to/file.itx"; loc = "SPD")

# Add any experiment-specific metadata needed by downstream analysis.
data = ARPESData(
    data;
    additional_metadata = Dict(
        :workfunction => 4.5,
    ),
)

# Adjust or reinterpret a physical axis without touching the intensity array.
data = convert_dim(data, :phi, x -> 0.5x)

# Convert the angular data to k-space.
kdata = k_conversion(data)
```

## Documentation map

```@contents
Pages = [
    "guide/data-model.md",
    "guide/io.md",
    "guide/analysis.md",
    "guide/k-space.md",
    "api.md",
]
Depth = 2
```

## Design overview

The source tree is intentionally layered:

- `src/ARPESData.jl` defines the main container and its metadata conventions.
- `src/io/` parses raw file formats and standardizes beamline-specific metadata.
- `src/dims.jl`, `src/transform.jl`, `src/filter.jl`, `src/derivative.jl` and `src/trapezoid.jl`
  provide dimension-aware processing utilities.
- `src/k_conv/` implements the k-space conversion pipeline.

If you are new to the package, start with
[Working with `ARPESData`](guide/data-model.md), then continue to
[Data loading](guide/io.md) and [k-space conversion](guide/k-space.md).
