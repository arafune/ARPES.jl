using DimensionalData
using DimensionalData.Lookups
using DimensionalData.Dimensions: label

export add_dim, shift_dim, negate_dim, rename_dim

"""
    add_dim(data::AbstractDimArray, dim_label, val::Real=NaN; unit::Union{Nothing, String}=nothing)

Add a new dimension to the `data` object with the specified `dim_label` and value `val`.
Optionally, a unit can be provided for the new dimension. Throws an error if the dimension already exists.

# Arguments
- `data::AbstractDimArray`: The AbstractDimArray object to modify.
- `dim_label::Symbol`: The label for the new dimension.
- `val::Real=NaN`: The value for the new dimension (default is NaN).
- `unit::Union{Nothing, String}=nothing`: Optional unit for the new dimension.

# Returns
- A new AbstractDimArray object with the added dimension.
"""
function add_dim(
    data::AbstractDimArray,
    dim_label::Symbol,
    val::Real;
    unit::Union{Nothing,String} = nothing,
)
    if hasdim(data, dim_label)
        throw(
            ArgumentError("A dimension with the name $dim_label already exists in $data."),
        )
    end
    new_data = reshape(parent(data), size(data)..., 1)
    newdim = Dim{dim_label}([val], metadata = Dict(:unit => unit))
    return rebuild(data; data = new_data, dims = (dims(data)..., newdim))
end

function add_dim(
    data::AbstractDimArray,
    dim_label::Symbol;
    unit::Union{Nothing,String} = nothing,
)
    if !haskey(metadata(data), dim_label)
        throw(
            ArgumentError(
                "A dimension with the name $dim_label does not exist in metadata of $data.",
            ),
        )
    end
    val = metadata(data)[dim_label]
    return add_dim(data, dim_label, val; unit = unit)
end

add_dim(
    data::AbstractDimArray,
    dim_label::String,
    val::Real;
    unit::Union{Nothing,String} = nothing,
) = add_dim(data, Symbol(dim_label), val; unit = unit)

function rename_dim(data, old::Symbol, new::Symbol)
    hasdim(data, old) || throw(ArgumentError("Dimension $old not found."))
    hasdim(data, new) && throw(ArgumentError("Dimension $new already exists."))
    ds = dims(data)
    i = findfirst(d -> name(d) == old, ds)

    d = ds[i]
    d2 = Dim{new}(lookup(d); metadata = metadata(d))

    return rebuild(data; dims = Base.setindex(ds, d2, i))
end

rename_dim(data, p::Pair{<:Symbol,<:Symbol}) = rename_dim(data, first(p), last(p))


"""
    shift_dim(dim::Dimension, shift_value)
    shift_dim(arpes_data::AbstractDimArray, dim_label::Symbol, shift_value)
    shift_dim(arpes_data::AbstractDimArray, dim::Type{<:DimensionalData.Dimension}, shift_value::Real)
    shift_dim(arpes_data::AbstractDimArray, dim::Dimension, shift_value)
    shift_dim(arpes_data::AbstractDimArray, dim_shift_pair...)

Shift the lookup values of a dimension or multiple dimensions in AbstractDimArray.

# Arguments
- `dim::Dimension`: The dimension to shift.
- `shift_value::Real`: The value to add to each element of the dimension's lookup.
- `arpes_data::AbstractDimArray`: The AbstractDimArray object containing dimensions.
- `dim_label::Symbol`: The name (type parameter) of the dimension to shift.
- `dim_shift_pair...`: A sequence of dimension labels (or Dimension objects) and shift values, e.g. `:phi, 1.0, :theta, 2.0`.

# Returns
- A new `Dimension` or `AbstractDimArray` with shifted lookup values.

# Examples
```julia
shift_dim(dim, 1.0)
shift_dim(arpes_data, :phi, 1.0)
shift_dim(arpes_data, :phi, 1.0, :theta, 2.0)
```
"""
function shift_dim(dim::DimensionalData.Dimension, shift_value::Real)
    new_lookup = _shift_lookup(lookup(dim), shift_value)
    return rebuild(dim; val = new_lookup)
end

function shift_dim(arpes_data::AbstractDimArray, dim_label::Symbol, shift_value::Real)
    idx = findfirst(d -> name(d) == dim_label, dims(arpes_data))

    if isnothing(idx)
        throw(
            ArgumentError("Dimension with name $dim_label not found in AbstractDimArray."),
        )
    end
    target_dim = dims(arpes_data)[idx]
    new_lookup = _shift_lookup(lookup(target_dim), shift_value)
    shifted_dim = rebuild(target_dim; val = new_lookup)
    new_dims = Base.setindex(dims(arpes_data), shifted_dim, idx)
    return rebuild(arpes_data; dims = new_dims)
end

function shift_dim(
    arpes_data::AbstractDimArray,
    dim::DimensionalData.Dimension,
    shift_value::Real,
)
    return shift_dim(arpes_data, name(dim), shift_value)
end

function shift_dim(
    arpes_data::AbstractDimArray,
    dim::Type{<:DimensionalData.Dimension},
    shift_value::Real,
)
    return shift_dim(arpes_data, name(dim), shift_value)
end

function shift_dim(arpes_data::AbstractDimArray, dim_shift_pair...)
    if isodd(length(dim_shift_pair))
        throw(
            ArgumentError(
                "dim_shift_pair must contain an even number of elements (dimension and shift pairs).",
            ),
        )
    end
    dims = dim_shift_pair[1:2:end]
    shifts = dim_shift_pair[2:2:end]
    for (dim, shift) in zip(dims, shifts)
        arpes_data = shift_dim(arpes_data, dim, shift)
    end
    return arpes_data
end

"""
  negate_dim(dim::Dimension)

Negate the lookup values of a dimension.
"""
function negate_dim(dim::DimensionalData.Dimension)
    new_lookup = _negate_lookup(lookup(dim))
    return rebuild(dim; val = new_lookup)
end

function change_energy_definition(arpes_data::ARPESData, new_definition::EnergyDefinition)
    #    energy_dim = findfirst(d -> name(d) == :eV, dims(arpes_data))
    #    if isnothing(energy_dim)
    #        throw(ArgumentError("Energy dimension with name :eV not found in AbstractDimArray."))
    #    end
    #    new_energy_dim = rebuild(
    #        dims(arpes_data)[energy_dim];
    #        val = lookup(dims(arpes_data)[energy_dim]),
    #        label = new_definition,
    #    )
    #    new_dims = map(
    #        (d, i) -> i == energy_dim ? new_energy_dim : d,
    #        dims(arpes_data),
    #        1:length(dims(arpes_data)),
    #    )
    #    return rebuild(arpes_data; dims = new_dims)
end

# -----------Private helper functions-----------

"""
    _shift_lookup(lookup, shift_value)

Return a new lookup array by adding `shift_value` to each element of `lookup`.

# Arguments
- `lookup`: The lookup array or range to shift.
- `shift_value`: The value to add to each element.

# Returns
- A shifted lookup array or range.
"""
function _shift_lookup(l::Lookup, shift)
    p = parent(l)
    return _shift_index(p, shift)
end

"""
    _negate_lookup(l::Lookup)

Negates the indices represented by the given `Lookup` object `l`.
Returns the negated index by calling `_negate_index` on the parent of `l`.
"""
function _negate_lookup(l::Lookup)
    p = parent(l)
    return _negate_index(p)
end

"""
    _process_index(index)

Process the given index for shifting. If the index is a range or array, return it as is.
Otherwise, throw an error.

# Arguments
- `index`: The index to process.

# Returns
- The processed index.

# Throws
- An error if the index type is not supported.
"""
function _shift_index(index::AbstractRange, shift)
    return range(index[1] + shift, step = step(index), length = length(index))
end

function _shift_index(index::AbstractArray, shift)
    return index .+ shift
end

"""
    _negate_index(index)

Negate the given index. If the index is a range or array, return the negated version.
"""
function _negate_index(index::AbstractRange)
    return range(- index[1], step = -step(index), length = length(index))
end

function _negate_index(index::AbstractArray)
    return -index
end
