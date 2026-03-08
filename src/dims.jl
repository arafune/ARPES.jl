using DimensionalData
using DimensionalData.Lookups

"""
    shift_dim(dim::Dimension, shift_value)

Return a new `Dimension` where the lookup values are shifted by `shift_value`.

# Arguments
- `dim::Dimension`: The dimension to shift.
- `shift_value`: The value to add to each element of the dimension's lookup.

# Returns
- A new `Dimension` with shifted lookup values.
"""
function shift_dim(dim::DimensionalData.Dimension, shift_value::Real)
    meta_data = metadata(dim)

    return rebuild(
        _shift_lookup(lookup(dim), shift_value);
        label = label(dim),
        metadata = meta_data,
    )
end

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
    new_index = _process_index(p, shift)
    return rebuild(l, new_index)
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
    return range(first(index) + shift, last(index) + shift, length = length(index))
end

function _process_index(index::AbstractArray, shift)
    return index .+ shift
end
