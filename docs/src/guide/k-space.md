```@meta
CurrentModule = ARPES
```

# k-space conversion

`k_conversion` transforms angle-domain ARPES data into momentum-space data.
This is one of the main analysis pipelines in the package.

## Expected input

The input must be `ARPESData` with:

- an `eV` dimension,
- a `phi` dimension,
- either a `psi` dimension or a dataset-level `:β` metadata entry, and
- metadata describing the experimental geometry.

The most important dataset-level metadata keys are:

| Key | Purpose |
| --- | --- |
| `:analyzer_configuration` | Selects the analyzer geometry mapping |
| `:workfunction` | Needed to convert energy to kinetic energy |
| `:hv` | Photon energy |
| `:β`, `:ξ`, `:δ`, `:χ` | Angular offsets and orientation parameters |

The `eV` dimension metadata must include `:energy_definition`.

Currently, the conversion supports:

- `BindingEnergy`
- `FinalStateEnergy`

## What the function does

Internally, the conversion pipeline:

1. validates the input metadata and dimensions,
2. converts the energy axis to the kinetic-energy convention needed for mapping,
3. interprets angular coordinates in radians,
4. determines or accepts target `kx` and `ky` ranges,
5. computes the inverse angle mapping for each target k-point, and
6. interpolates the intensity onto the momentum grid.

The result is returned as a new `ARPESData` with momentum-space dimensions.

## Output

The output keeps the `ARPESData` abstraction and carries over the dataset metadata.
The main axes are:

- `kx`
- `ky`
- `eV`

Optional `kz` support is part of the API surface, but the current implementation is mainly
centered on the `kx`/`ky` workflow.

## Minimal example

```julia
using ARPES

data = load(
    "path/to/file.itx";
    extra_metadata = Dict(
        :workfunction => 4.5,
        :ξ => 0.0,
        :δ => 0.0,
        :χ => 0.0,
    ),
)

kdata = k_conversion(data)
```

If you need a specific output grid, pass the target ranges explicitly:

```julia
kdata = k_conversion(
    data;
    kx_range = range(-1.0, 1.0; length = 241),
    ky_range = range(-1.0, 1.0; length = 241),
)
```

## Offsets and conventions

The keywords `β0`, `ξ0`, `δ0`, and `χ0` provide additional offsets on top of the metadata
stored in the dataset.

This is useful when:

- refining alignment after loading,
- comparing scans with slightly different mounting conditions, or
- exploring sensitivity to angular calibration.

## Practical advice

Before calling `k_conversion`, check the following:

- the `eV` dimension has the right energy definition,
- the dataset contains a valid analyzer configuration,
- `:workfunction` is present, and
- the angular metadata matches the instrument convention of the imported data.

In most workflows, the most reliable pattern is:

1. load and standardize the data,
2. inspect or adjust metadata,
3. optionally smooth or align the data, and
4. run `k_conversion`.

See the [API reference](../api.md) for the `k_conversion` signature and docstring.
