module Location
using ..IO: LocationLoader, SPDLoader, GenericLoader
export canonical_dim, dim_alias, default_dim_map

dim_alias(::LocationLoader) = Dict{Function,Vector{Symbol}}()
default_dim_map(::LocationLoader) = Dict{Symbol,Function}()

function canonical_dim(loader::LocationLoader, name::Symbol)
    aliases = dim_alias(loader)

    for (ctor, aliases) in aliases
        if name in aliases
            return ctor
        end
    end
    return get(default_dim_map(loader), name, nothing)
end

include("spd.jl")
export load_data
end
