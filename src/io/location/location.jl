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
    angle_Shin_convention(::Type{<:LocationLoader}) -> Dict
    angle_Shin_convention(loader::LocationLoader) -> Dict

Returns the angle standardization configuration for a given `LocationLoader` instance by dispatching to its type.

# Arguments
- `loader::Type{<:LocationLoader}`: The loader type.
- `loader::LocationLoader`: The loader instance.

# Returns
- `Dict{Symbol,Union{Dict{Symbol,Function},Number,Vector{Dict{Symbol,Function}}}}`: An empty dictionary.
"""
angle_Shin_convention(::Type{<:LocationLoader}) =
    Dict{Symbol,Union{Dict{Symbol,Function},Number,Vector{Dict{Symbol,Function}}}}()
angle_Shin_convention(loader::LocationLoader) = angle_Shin_convention(typeof(loader))



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

"""
  convert_Shin_convention(metadata_dict::Dict, conversion_rule::Dict) -> Dict{Symbol,Any}

Convert metadata keys according to the Shin convention.

Given a `metadata_dict` containing raw metadata and a `conversion_rule` dict specifying how to map and resolve keys,
this function applies the conversion rules and returns a new dictionary with keys as specified in `conversion_rule`
and values resolved using `_resolve_Shin_rule`. Only keys with non-`nothing` resolved values are included in the result.

# Arguments
- `metadata_dict::Dict`: The input metadata dictionary.
- `conversion_rule::Dict`: A dictionary mapping target keys to conversion rules.

# Returns
- `Dict{Symbol,Any}`: A dictionary with converted keys and resolved values.
"""
function convert_Shin_convention(metadata_dict::Dict, conversion_rule::Dict)

    result = Dict{Symbol,Any}()
    for (target_key, rule) in conversion_rule
        value = _resolve_Shin_rule(metadata_dict, target_key, rule)
        value !== nothing && (result[target_key] = value)
    end

    return result
end

"""
 _resolve_Shin_rule: Resolves a rule against metadata.

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
function _resolve_Shin_rule(metadata::Dict, _::Symbol, rule::Symbol)
    return get(metadata, rule, nothing)
end

function _resolve_Shin_rule(metadata::Dict, _::Symbol, rule::Dict)
    for (source_key, f) in rule
        if haskey(metadata, source_key)
            return f(metadata[source_key])
        end
    end
    return nothing
end

function _resolve_Shin_rule(metadata::Dict, target_key::Symbol, rule::Function)
    haskey(metadata, target_key) || return nothing
    return rule(metadata[target_key])
end

function _resolve_Shin_rule(metadata::Dict, target_key::Symbol, rule::Vector)
    for subrule in rule
        value = _resolve_Shin_rule(metadata, target_key, subrule)
        value !== nothing && return value
    end
    return nothing
end

function _resolve_Shin_rule(_::Dict, _::Symbol, rule)
    return rule
end

include("spd.jl")
export load_data
end
