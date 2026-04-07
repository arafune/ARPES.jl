"""
  Location

A module for handling location-specific data loading and standardization.
This module provides functionality for converting ARPESData from DimArray
to standardized form according to canonical dimension mappings and metadata conversion rules.

The `Location` module includes functions for:

- Mapping dimension constructor functions to their alias symbols for different `LocationLoader` types.
- Mapping dimension alias symbols to their constructor functions for different `LocationLoader` types.
- Standardizing dimensions of raw `DimArray` data according to canonical dimension mappings.
- Converting metadata keys according to specified conversion rules.

The module also includes the `spd.jl` file, which contains the implementation for loading SPD data.
Additional location-specific loaders can be implemented and included in this module as needed.


In the file of Location module, `load_data` function must be fefined for each specific loader type, 
in which the function defined in Format module should be called to read the data in its original form,
and then the `to_standardize` function should be used to convert it to the standardized form.
Further, The `canonical_dim`, `dim_alias`, `default_dim_map`, and `metadata_convert_rule` functions should be
implemented to provide the necessary mappings and conversion rules for that loader.
"""

module Location
using DimensionalData
using DimensionalData.Dimensions.Lookups: NoMetadata
using ARPES.IO: LocationLoader, SPDLoader, GenericLoader
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
    metadata_convert_rule(::Type{<:LocationLoader}) -> Dict
    metadata_convert_rule(loader::LocationLoader) -> Dict

Returns the angle standardization configuration for a given `LocationLoader` instance by dispatching to its type.

# Arguments
- `loader::Type{<:LocationLoader}`: The loader type.
- `loader::LocationLoader`: The loader instance.

# Returns
- `Dict{Symbol,Union{Dict{Symbol,Function},Number,Vector{Dict{Symbol,Function}}}}`: An empty dictionary.
"""
metadata_convert_rule(::Type{<:LocationLoader}) =
    Dict{Symbol,Union{Dict{Symbol,Function},Number,Vector{Dict{Symbol,Function}}}}()
metadata_convert_rule(loader::LocationLoader) = metadata_convert_rule(typeof(loader))



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

    metadata_dict = Dict(metadata(raw))
    converted_metadata = standardize_metadata(metadata_dict, metadata_convert_rule(loader))
    new_metadata = merge(metadata_dict, converted_metadata)

    return rebuild(raw; dims = new_dims, metadata = new_metadata)
end

"""
Finds the canonical dimension constructor for a given alias symbol.

Searches through the dimension aliases and returns the constructor function if the alias is found.
If not found, returns the constructor from the default dimension map or `nothing` if not present.

# Arguments
- `loader::Type{<:LocationLoader}`: The loader type.
- `name::Symbol`: The alias symbol to resolve.

# Returns
- `Symbol` or `Nothing`: The canonical name of angle defined by RevSciInstruments or `nothing` if not found.
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

"""
  standardize_metadata(metadata_dict::Dict, conversion_rule::Dict) -> Dict{Symbol,Any}

Convert metadata keys according to conversion_rule convention.

Given a `metadata_dict` containing raw metadata and a `conversion_rule` dict specifying how to map and resolve keys,
this function applies the conversion rules and returns a new dictionary with keys as specified in `conversion_rule`
and values resolved using `_normalize_metadata`.
Only keys with non-`nothing` resolved values are included in the result.

# Arguments
- `metadata_dict::Dict`: The input metadata dictionary.
- `conversion_rule::Dict`: A dictionary mapping target keys to conversion rules.

# Returns
- `Dict{Symbol,Any}`: A dictionary with converted keys and resolved values.
"""
function standardize_metadata(metadata_dict::Dict, conversion_rule::Dict)

    result = Dict{Symbol,Any}()
    for (target_key, rule) in conversion_rule
        value = _normalize_metadata(metadata_dict, target_key, rule)
        value !== nothing && (result[target_key] = value)
    end

    return result
end

"""
 _normalize_metadata: Resolves a rule against metadata.

 This function is overloaded to handle different types of `rule`:
   1. Symbol: Returns the value for the symbol key in metadata, or `nothing` if not present.
   2. Dict: Iterates over key-function pairs, applies the function to the value if the key exists.
   3. Function: Applies the function to the value of the target key if it exists, otherwise returns `nothing`.
   4. Vector: Tries each subrule in order, returning the first non-nothing result.
   5. Value: Returns the rule itself.

 Arguments:
   metadata::Dict : Dictionary containing metadata.
   target_key::Symbol : The key for which the rule is being resolved (used for function rules).
   rule           : Rule to resolve (Symbol, Dict, Vector, or any value).

 Returns:
   The resolved value according to the rule type, or `nothing` if not found.
"""
function _normalize_metadata(metadata::Dict, _::Symbol, rule::Symbol)
    return get(metadata, rule, nothing)
end

function _normalize_metadata(metadata::Dict, _::Symbol, rule::Dict)
    for (source_key, f) in rule
        if haskey(metadata, source_key)
            return f(metadata[source_key])
        end
    end
    return nothing
end

function _normalize_metadata(metadata::Dict, target_key::Symbol, rule::Function)
    haskey(metadata, target_key) || return nothing
    return rule(metadata[target_key])
end

function _normalize_metadata(metadata::Dict, target_key::Symbol, rule::Vector)
    for subrule in rule
        value = _normalize_metadata(metadata, target_key, subrule)
        value !== nothing && return value
    end
    return nothing
end

function _normalize_metadata(_::Dict, _::Symbol, rule)
    return rule
end

include("spd.jl")
export load_data
end
