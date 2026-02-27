module Location
using DimensionalData
using ..IO: LocationLoader, SPDLoader, GenericLoader
export canonical_dim, dim_alias, default_dim_map

"""
Returns a dictionary mapping dimension constructor functions to their alias symbols for a LocationLoader.
"""
dim_alias(::Type{<:LocationLoader}) = Dict{Function,Vector{Symbol}}()
dim_alias(loader::LocationLoader) = dim_alias(typeof(loader))

"""
Returns a dictionary mapping dimension alias symbols to their constructor functions for a LocationLoader.
"""
default_dim_map(::Type{<:LocationLoader}) = Dict{Symbol,Function}()
default_dim_map(loader::LocationLoader) = default_dim_map(typeof(loader))

"""
    angle_standardization(::Type{<:LocationLoader}) -> Dict

Returns an empty dictionary for angle standardization configuration for a given `LocationLoader` type.
# Returns
- `Dict{Symbol,Union{Dict{Symbol,Function},Number,Vector{Dict{Symbol,Function}}}}`: An empty dictionary.
"""
angle_standardization(::Type{<:LocationLoader}) =
    Dict{Symbol,Union{Dict{Symbol,Function},Number,Vector{Dict{Symbol,Function}}}}()
"""
    angle_standardization(loader::LocationLoader) -> Dict

Returns the angle standardization configuration for a given `LocationLoader` instance by dispatching to its type.

# Arguments
- `loader::LocationLoader`: The loader instance.

# Returns
- `Dict{Symbol,Union{Dict{Symbol,Function},Number,Vector{Dict{Symbol,Function}}}}`: The configuration dictionary.
"""
angle_standardization(loader::LocationLoader) = angle_standardization(typeof(loader))



"""
    to_standardize(loader::Type{<:LocationLoader}, raw::DimArray)

Standardizes the dimensions of a raw `DimArray` according to the canonical dimension mapping
for the given `LocationLoader` type.

Each dimension in `raw` is checked against the canonical mapping. If a canonical constructor
is found for the dimension name, it is used to reconstruct the dimension; otherwise, the original
dimension is retained.

# Arguments
- `loader::Type{<:LocationLoader}`: The loader type for which to standardize dimensions.
- `raw::DimArray`: The raw DimArray to be standardized.

# Returns
A new `DimArray` with standardized dimensions.
"""
function to_standardize(loader::Type{<:LocationLoader}, raw::DimArray)
    new_dims = map(dims(raw)) do d
        raw_name = name(d)
        ctor=canonical_dim(loader, raw_name)
        isnothing(ctor) ? d : ctor(val(d))
    end


    return rebuild(raw; dims = new_dims)
end

"""
Finds the canonical dimension constructor for a given alias symbol.

Searches through the dimension aliases and returns the constructor function if the alias is found.
If not found, returns the constructor from the default dimension map or `nothing` if not present.

# Arguments
- `loader::Type{<:LocationLoader}`: The loader type.
- `name::Symbol`: The alias symbol to resolve.

# Returns
- `Function` or `Nothing`: The canonical constructor function or `nothing` if not found.
"""
function canonical_dim(loader::Type{<:LocationLoader}, name::Symbol)
    alias_map = dim_alias(loader)

    for (ctor, alias_list) in alias_map
        if name in alias_list
            return ctor
        end
    end

    return get(default_dim_map(loader), name, nothing)
end

canonical_dim(loader::LocationLoader, name::Symbol) = canonical_dim(typeof(loader), name)


function angles_for_k_conversion(loaader::Type{<:LocationLoader}, metadata) end



include("spd.jl")
export load_data
end
