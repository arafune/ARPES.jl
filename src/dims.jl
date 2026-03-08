using DimensionalData
using DimensionalData.Lookups
using DimensionalData.Dimensions: label

export shift_dim

"""
    shift_dim(dim::Dimension, shift_value)
    shift_dim(arpes_data::ARPESData, dim_label::Symbol, shift_value)
    shift_dim(arpes_data::ARPESData, dim::Type{<:DimensionalData.Dimension}, shift_value::Real)
    shift_dim(arpes_data::ARPESData, dim::Dimension, shift_value)
    shift_dim(arpes_data::ARPESData, dim_shift_pair...)

Shift the lookup values of a dimension or multiple dimensions in ARPESData.

# Arguments
- `dim::Dimension`: The dimension to shift.
- `shift_value::Real`: The value to add to each element of the dimension's lookup.
- `arpes_data::ARPESData`: The ARPESData object containing dimensions.
- `dim_label::Symbol`: The name (type parameter) of the dimension to shift.
- `dim_shift_pair...`: A sequence of dimension labels (or Dimension objects) and shift values, e.g. `:phi, 1.0, :theta, 2.0`.

# Returns
- A new `Dimension` or `ARPESData` with shifted lookup values.

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

function shift_dim(arpes_data::ARPESData, dim_label::Symbol, shift_value::Real)
    idx = findfirst(d -> name(d) == dim_label, dims(arpes_data))

    if isnothing(idx)
        throw(ArgumentError("Dimension with name $dim_label not found in ARPESData."))
    end
    target_dim = dims(arpes_data)[idx]
    new_lookup = _shift_lookup(lookup(target_dim), shift_value)
    shifted_dim = rebuild(target_dim; val = new_lookup)
    new_dims = Base.setindex(dims(arpes_data), shifted_dim, idx)
    return rebuild(arpes_data; dims = new_dims)
end

function shift_dim(arpes_data::ARPESData, dim::DimensionalData.Dimension, shift_value::Real)
    return shift_dim(arpes_data, name(dim), shift_value)
end

function shift_dim(
    arpes_data::ARPESData,
    dim::Type{<:DimensionalData.Dimension},
    shift_value::Real,
)
    return shift_dim(arpes_data, name(dim), shift_value)
end

function shift_dim(arpes_data::ARPESData, dim_shift_pair...)
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

function change_energy_definition(arpes_data::ARPESData, new_definition::EnergyDefinition)
    #    energy_dim = findfirst(d -> name(d) == :eV, dims(arpes_data))
    #    if isnothing(energy_dim)
    #        throw(ArgumentError("Energy dimension with name :eV not found in ARPESData."))
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
    return _process_index(p, shift)
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
function _process_index(index::AbstractRange, shift)
    return range(index[1] + shift, step = step(index), length = length(index))
end

function _process_index(index::AbstractArray, shift)
    return index .+ shift
end
