# Copilot instructions for ARPES.jl

## Development commands

### Main package environment

- Instantiate and build the package environment:

  ```bash
  julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
  ```

- Run the full test suite:

  ```bash
  julia --project=. -e 'using Pkg; Pkg.test()'
  ```

- Run a single test file from the repository root:

  ```bash
  julia --project=. -e 'include("test/io/io.jl")'
  ```

  This pattern also works for other test files under `test/` because the suite is organized as plain Julia files that can be included directly.

### Documentation environment

- Build the docs with the dedicated docs project:

  ```bash
  julia --project=docs docs/make.jl
  ```

- If the docs environment has not been instantiated yet:

  ```bash
  julia --project=docs -e 'using Pkg; Pkg.instantiate()'
  ```

### Formatting

- The repository uses JuliaFormatter's `blue` style (`.JuliaFormatter.toml`):

  ```bash
  julia --project=. -e 'using JuliaFormatter; format("src")'
  ```

## High-level architecture

- `src/ARPES.jl` is the wiring point for the package. It includes the core type definitions, array/data helpers, filter operations, I/O stack, dimension helpers, and k-space conversion code, then re-exports the main public API (`ARPESData`, dimension labels, `load`, and `k_conversion`).

- `ARPESData` is the central abstraction. It wraps a `DimensionalData.AbstractDimArray`, preserves dimension-aware indexing/metadata behavior, and delegates the standard DimensionalData interface through methods like `parent`, `dims`, `metadata`, and `rebuild`. Most higher-level functionality assumes data stays in this dimension-aware form rather than being converted to plain arrays.

- The I/O pipeline is intentionally layered:
  - `src/io/formats/` parses raw file formats into a mostly faithful `DimArray`.
  - `src/io/location/` standardizes beamline/location-specific dimension names and metadata into the package's canonical ARPES conventions.
  - `src/io/registry.jl` selects a loader either from `loc=` or by sniffing file content / extension.
  - The standardized result is finally wrapped as `ARPESData`.

- `SPDLoader` is the concrete example of that pipeline today. It reads raw ITX data, maps aliases like `:theta_y` or `:energy_channel` onto canonical dimensions such as `phi`, `eV`, and `detector_ch`, normalizes metadata, and then constructs `ARPESData`. `.sp2` loader selection exists in the registry, but parsing is not implemented yet.

- `src/k_conv/` assumes the data is already standardized. `k_conversion` reads canonical dimensions and metadata from `ARPESData`, performs preprocessing + interpolation + postprocessing, and returns a new `ARPESData` on momentum-space axes (`kx`, `ky`, `kz` as available).

- `test/runtests.jl` is the suite entrypoint and mirrors the source layout: top-level feature groups include smaller files from `test/`, `test/io/`, and `test/k_conv/`. Tests rely on fixture files in `testdata/`, which is a git submodule.

## Repository-specific conventions

- Preserve the `DimensionalData` model. When changing dimensions or metadata, follow the existing pattern of using `dims`, `lookup`, `metadata`, `parent`, and `rebuild` instead of manually unpacking and rebuilding arrays ad hoc.

- Treat metadata as part of the type contract, not incidental data. In particular:
  - dataset-level metadata commonly carries keys such as `:intensity_unit`, `:analyzer_configuration`, `:hv`, `:workfunction`, and angular offsets like `:β`, `:ξ`, `:δ`, `:χ`
  - the `eV` dimension metadata is expected to include `:energy_definition`
  - `k_conversion` depends on those canonical metadata keys being present and meaningful

- Keep raw-format parsing separate from location-specific standardization. New loaders should follow the same split used by `Format` and `Location`: parse the file into a raw `DimArray`, then standardize dimension aliases and metadata through `Location` methods before constructing `ARPESData`.

- If you add a new location loader, wire all of the following surfaces together:
  - define a `LocationLoader` subtype in `src/io/types.jl`
  - register it in `LOCATION_REGISTRY`
  - implement `dim_alias`, `default_dim_map`, `metadata_convert_rule`, and `load_data`
  - add focused tests under `test/io/`

- The test suite uses reusable helpers from `test/fixture.jl`. Test files often guard fixture includes with `if !@isdefined(...)` so they can be run both from `test/runtests.jl` and directly as single included files. Preserve that pattern when adding targeted tests.

- `test/` generally mirrors `src/` by feature area (`io`, `k_conv`, core array helpers). Follow that structure instead of creating unrelated one-off test locations, and remember to include any new test file from `test/runtests.jl`.

- Keep repository-facing docstrings and code comments in English. The docs build is minimal and based on `Documenter` autodocs, so public API docstrings feed directly into generated documentation.
