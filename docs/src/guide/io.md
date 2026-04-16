```@meta
CurrentModule = ARPES
```

# Data loading and standardization

`ARPES.jl` separates raw file parsing from location-specific standardization.
This keeps file-format code focused on reading data, while beamline-specific code maps
dimensions and metadata into the package's canonical ARPES model.

## The `load` entry point

The main user-facing entry point is `load`:

```julia
using ARPES

data = load("path/to/file.itx")
```

You can optionally force a location-specific loader:

```julia
data = load("path/to/file.itx"; loc = "SPD")
```

You can also inject extra metadata during loading:

```julia
data = load(
    "path/to/file.itx";
    extra_metadata = Dict(:workfunction => 4.5, :sample_temperature => 20.0),
)
```

This is especially useful when a raw file does not contain all metadata required for
downstream analysis.

## How loader selection works

The loading pipeline is:

1. `load` calls `select_loader`.
2. The registry either uses the explicit `loc` argument or detects a loader from the file.
3. The selected loader reads the raw file format.
4. Location-specific standardization renames dimensions and converts metadata.
5. The standardized result is returned as `ARPESData`.

## Supported formats and locations

The current source tree includes:

- ITX parsing through `read_itx`,
- loader detection for `.sp2` files (parsing not yet implemented), and
- an `SPDLoader` that standardizes dimensions and metadata into canonical ARPES names.

In other words, SPD is the current concrete location-specific loader in the repository, and
it is the main path exercised by the user-facing examples.

When no explicit location is supplied, loader detection uses filename extensions and, for
ITX files, header inspection.

## SPD standardization

The SPD pipeline is the concrete reference implementation in the repository.
It standardizes:

- angle-like raw dimensions to `phi`,
- energy-like raw dimensions to `eV`,
- detector-channel dimensions to `detector_ch`,
- metadata such as photon energy to canonical keys like `:hv`.

It also derives angular metadata used later by `k_conversion`.

For SPD data, the loader sets:

- `:analyzer_configuration` to `TypeI`,
- the energy definition to `FinalStateEnergy`, and
- the intensity unit to `Counts` or `CPS` based on the raw metadata.

## When to add metadata manually

Not every file carries the full experimental context needed for advanced analysis.
In particular, k-space conversion requires metadata such as:

- `:workfunction`
- `:hv`
- `:β`, `:ξ`, `:δ`, `:χ` when they are not represented directly by dimensions

If any of these are absent after loading, add them explicitly before conversion.

## Extending the I/O layer

The source code is designed so new location loaders follow the same pattern:

1. parse a raw file into a dimensional array,
2. define dimension aliases and default mappings,
3. define metadata conversion rules, and
4. wrap the standardized result as `ARPESData`.

This split makes it easier to support new beamlines without changing the rest of the
analysis pipeline.

See the [API reference](../api.md) for the `load` signature and docstring.
