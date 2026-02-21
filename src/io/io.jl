module IO
include("types.jl")
include("registry.jl")
include("location/location.jl")
include("formats/format.jl")

function load(fpath::String; loc = nothing)
    #  resolved_path = resolve_path(fpath)

    loader_type = select_loader(resolved_path, loc)

    return load_data(loader_type, path)
end
end

