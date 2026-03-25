```@meta
CurrentModule = ARPES
```

# Working with `ARPESData`

`ARPESData` is the central abstraction in `ARPES.jl`.
It wraps a `DimensionalData.AbstractDimArray`, so the array values, axes, and metadata
stay coupled during indexing and reconstruction.

## Core idea

An `ARPESData` object stores:

- an N-dimensional intensity array,
- named dimensions such as `eV`, `phi`, `psi`, or `kx`,
- an optional dataset name, and
- metadata required for analysis pipelines.

This lets the package operate on physical axes instead of anonymous array indices.

## Canonical dimensions

The package defines canonical dimension labels for common ARPES workflows:

| Dimension | Meaning |
| --- | --- |
| `eV` | Energy axis |
| `phi` | Emission angle parallel to the analyzer slit |
| `psi` | Emission angle perpendicular to the analyzer slit |
| `kx`, `ky`, `kz` | Momentum-space axes in reciprocal space |
| `detector_ch`, `ch2` | Detector-channel axes |
| `delay` | Pump-probe delay axis |
| `spin` | Spin-resolved axis |

In practice, the same dataset may move between angle-space dimensions (`phi`, `psi`) and
momentum-space dimensions (`kx`, `ky`, `kz`) while preserving the same metadata model.

## Metadata contract

Metadata is not treated as optional decoration.
Several operations, especially `k_conversion`, expect canonical keys to be present and
meaningful.

Common dataset-level metadata keys are:

| Key | Meaning |
| --- | --- |
| `:intensity_unit` | Intensity unit such as `Counts` or `CPS` |
| `:analyzer_configuration` | Analyzer geometry such as `TypeI` or `TypeII` |
| `:hv` | Photon energy in eV |
| `:workfunction` | Sample work function in eV |
| `:β`, `:ξ`, `:δ`, `:χ` | Angular offsets used by k-space conversion |

The `eV` dimension should also carry dimension metadata with the key
`:energy_definition`, whose value is one of:

- `BindingEnergy`
- `FinalStateEnergy`
- `KineticEnergy`
- `IntermediateEnergy`

At present, `k_conversion` supports `BindingEnergy` and `FinalStateEnergy`.

## Constructing data

There are two common construction paths.

### From raw arrays and dimensions

```julia
using ARPES, DimensionalData

energy = eV(range(-0.2, 0.2; length = 200), metadata = Dict(:energy_definition => BindingEnergy))
angle = phi(range(-12.0, 12.0; length = 181))

raw = rand(length(energy), length(angle))

data = ARPESData(
    raw,
    (energy, angle);
    metadata = Dict(
        :intensity_unit => Counts,
        :analyzer_configuration => TypeI,
        :hv => 21.2,
        :workfunction => 4.5,
        :β => 0.0,
        :ξ => 0.0,
        :δ => 0.0,
        :χ => 0.0,
    ),
)
```

### From an existing dimensional array

```julia
dimarray = DimArray(raw, (energy, angle))

data = ARPESData(
    dimarray;
    intensity_unit = Counts,
    analyzer_config = TypeI,
    energy_def = BindingEnergy,
    additional_metadata = Dict(:hv => 21.2, :workfunction => 4.5),
)
```

This constructor injects `:intensity_unit` and `:analyzer_configuration` at the dataset
level, and writes `:energy_definition` into the `eV` dimension metadata.

## Dimension-aware editing

Several utilities modify dimensions while preserving the dimensional-array model:

- `add_dim` adds a scalar dimension,
- `rename_dim` renames an axis,
- `shift_dim` shifts lookup values without changing array contents, and
- `negate_dim` flips the sign of a dimension lookup.

These are useful when calibrating offsets, aligning scans, or adapting imported data to
canonical conventions.

## Why this matters

Because the package keeps dimensions and metadata attached to the data:

- transforms can preserve the physical meaning of each axis,
- file loaders can standardize different beamline conventions into one schema, and
- `k_conversion` can use the same metadata contract regardless of input format.

See the [API reference](../api.md) for the constructor docstring and signature details.
