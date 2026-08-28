using DimensionalData
using DimensionalData.Lookups
using DimensionalData.Dimensions: label

export add_dim, convert_dim, shift_dim, negate_dim, rename_dim

"""
    add_dim(data::AbstractDimArray, dim_label, val::Real=NaN; unit::Union{Nothing, String, Symbol}=nothing)

Add a new dimension to the `data` object with the specified `dim_label` and value `val`.
Optionally, a unit can be provided for the new dimension. Throws an error if the dimension already exists.

# Arguments
- `data::AbstractDimArray`: The AbstractDimArray object to modify.
- `dim::Symbol`: The label for the new dimension.
- `val::Real=NaN`: The value for the new dimension (default is NaN).
- `unit::Union{Nothing, String, Symbol}=nothing`: Optional unit for the new dimension.

# Returns
- A new AbstractDimArray object with the added dimension.
"""
function add_dim(
    data::AbstractDimArray,
    dim::Symbol,
    val::Real;
    unit::Union{Nothing,String,Symbol} = nothing,
)
    if hasdim(data, dim)
        throw(ArgumentError("A dimension with the name $dim already exists in $data."))
    end
    new_data = reshape(parent(data), size(data)..., 1)
    newdim = Dim{dim}([val], metadata = Dict(:unit => unit))
    return rebuild(data; data = new_data, dims = (dims(data)..., newdim))
end

function add_dim(
    data::AbstractDimArray,
    dim::Symbol;
    unit::Union{Nothing,String,Symbol} = nothing,
)
    if !haskey(metadata(data), dim)
        throw(
            ArgumentError(
                "A dimension with the name $dim does not exist in metadata of $data.",
            ),
        )
    end
    val = metadata(data)[dim]
    return add_dim(data, dim, val; unit = unit)
end

add_dim(
    data::AbstractDimArray,
    dim::String,
    val::Real;
    unit::Union{Nothing,String,Symbol} = nothing,
) = add_dim(data, Symbol(dim), val; unit = unit)


function add_dim(
    data::AbstractDimArray,
    dim::DimensionalData.Dimension,
    val::Real;
    unit::Union{Nothing,String,Symbol} = nothing,
)
    dim_label = DimensionalData.name(dim)
    return add_dim(data, dim_label, val; unit = unit)
end

function add_dim(
    data::AbstractDimArray,
    dim::DimensionalData.Dimension;
    unit::Union{Nothing,String,Symbol} = nothing,
)
    dim_label = DimensionalData.name(dim)

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

"""
    convert_dim(dim::DimensionalData.Dimension, f::Function)
    convert_dim(data::AbstractDimArray, dim::Union{DimensionalData.Dimension,Symbol}, f::Function)

Apply a coordinate transform `f` to the lookup values of a dimension and return a new
dimension with preserved metadata. When called on an `AbstractDimArray`, only the selected
axis lookup is transformed and the array is rebuilt with the same parent data and metadata.

# Arguments
- `dim::DimensionalData.Dimension`: An existing dimension object whose lookup values are
  transformed. For example, this can be a value such as `Dim{:x}(0:0.5:10)` taken from a
  dimensional array or constructed directly.
- `data::AbstractDimArray`: A dimensional array whose selected dimension should be
  transformed in place via rebuilding.
- `dim::Union{DimensionalData.Dimension,Symbol}`: The dimension to transform in `data`,
  identified either by the dimension object itself or by its name.
- `f::Function`: A function applied elementwise to the lookup values of `dim`.

# Returns
- A new `DimensionalData.Dimension` with transformed coordinates and the same metadata as
  `dim`. The original dimension name is preserved.

If the transformed coordinates remain equally spaced, as determined by the shared
`_is_equal_spacing` helper from `utils.jl`, the returned dimension uses a range lookup so
`step` is still available. Otherwise the transformed coordinates are stored as a vector
lookup.

# Examples
```julia
dim = Dim{:x}(0:0.5:10)

convert_dim(dim, x -> 2x + 1)
```
"""
function convert_dim(dim::DimensionalData.Dimension, f::Function)
    converted_values = f.(parent(lookup(dim)))
    new_lookup = _converted_lookup(converted_values)
    return Dim{name(dim)}(new_lookup; metadata = metadata(dim))
end

function convert_dim(
    data::AbstractDimArray,
    dim::Union{DimensionalData.Dimension,Symbol},
    f::Function,
)
    if !hasdim(data, dim)
        throw(ArgumentError("Dimension with name $dim not found in AbstractDimArray."))
    end
    dim = dims(data, dim)
    idx = dimnum(data, dim)

    converted = convert_dim(dim, f)
    new_dims = Base.setindex(dims(data), converted, idx)
    return rebuild(data; dims = new_dims)
end

function convert_dim(_::DimensionalData.Dimension, f)
    throw(ArgumentError("f must be a function."))
end

function convert_dim(_::AbstractDimArray, _::Union{DimensionalData.Dimension,Symbol}, f)
    throw(ArgumentError("f must be a function."))
end

"""
    rename_dim(data::AbstractDimArray, old::Symbol, new::Symbol)
    rename_dim(data::AbstractDimArray, p::Pair{<:Symbol,<:Symbol})

Rename a dimension in a dimensional array while preserving its lookup values and metadata.

# Arguments
- `data::AbstractDimArray`: Input dimensional array.
- `old::Symbol`: Existing dimension label.
- `new::Symbol`: New dimension label.
- `p::Pair`: Convenience form such as `:theta => :phi`.

# Returns
- A rebuilt dimensional array with the renamed dimension.

# Errors
- Throws `ArgumentError` if `old` does not exist or if `new` already exists.
"""
function rename_dim(data, old::Symbol, new::Symbol)
    hasdim(data, old) || throw(ArgumentError("Dimension $old not found."))
    hasdim(data, new) && throw(ArgumentError("Dimension $new already exists."))
    ds = dims(data)
    i = dimnum(ds, old)

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
    if !hasdim(arpes_data, dim_label)
        throw(
            ArgumentError("Dimension with name $dim_label not found in AbstractDimArray."),
        )
    end

    idx = dimnum(arpes_data, dim_label)

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
    dim_labels = dim_shift_pair[1:2:end]
    shifts = dim_shift_pair[2:2:end]
    for (dim, shift) in zip(dim_labels, shifts)
        arpes_data = shift_dim(arpes_data, dim, shift)
    end
    return arpes_data
end

"""
    negate_dim(dim::Dimension)
    negate_dim(data::AbstractDimArray, dim_label::Symbol)
    negate_dim(data::AbstractDimArray, dim::DimensionalData.Dimension)

Negate the lookup values of a dimension.
"""
function negate_dim(dim::DimensionalData.Dimension)
    old_lookup = lookup(dim)
    new_lookup = _negate_lookup(lookup(dim))

    @debug "negate_dim_check" length(parent(old_lookup)) length(new_lookup)

    negated_dim = rebuild(dim; val = new_lookup)

    @debug "negate_dim after rebuild" length(parent(lookup(negated_dim)))
    return negated_dim
end

function negate_dim(data::AbstractDimArray, dim_label::Symbol)
    if !hasdim(data, dim_label)
        throw(
            ArgumentError("Dimension with name $dim_label not found in AbstractDimArray."),
        )
    end

    idx = dimnum(data, dim_label)
    negated = negate_dim(dims(data)[idx])
    return rebuild(data; dims = Base.setindex(dims(data), negated, idx))
end

function negate_dim(data::AbstractDimArray, dim::DimensionalData.Dimension)
    return negate_dim(data, name(dim))
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
    _shift_index(index, shift)

Return a new index by adding `shift` to each element of `index`.
"""
function _shift_index(index::AbstractRange, shift)
    return range(index[1] + shift, step = step(index), length = length(index))
end

function _shift_index(index::AbstractArray, shift)
    return index .+ shift
end

function _converted_lookup(values::AbstractVector)
    if length(values) ≥ 2 && _is_equal_spacing(values)
        spacing = (values[end] - values[1]) / (length(values) - 1)
        return range(values[1], step = spacing, length = length(values))
    end
    return collect(values)
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
