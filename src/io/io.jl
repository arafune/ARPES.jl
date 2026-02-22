"""
    IO

Module for input/output operations in ARPES.jl.
Includes functions and types for loading and saving data,
handling file formats, and managing file locations.
"""
module IO
include("types.jl")
include("registry.jl")
include("location/location.jl")
include("formats/format.jl")

"""
    load(fpath::String; loc = Union{String, nothing})

Load data from the given file path.
- `fpath`: Path to the file.
- `loc`: Optional location argument (String or nothing).

Returns parsed data.

# Example
```julia
data = IO.load("data/file.itx", loc="spd")
```

"""
function load(fpath::String; loc = Union{String,nothing})
    #  resolved_path = resolve_path(fpath)

    loader_type = select_loader(resolved_path, loc)

    return load_data(loader_type, path)
end
export load
end

