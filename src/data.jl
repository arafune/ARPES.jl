using Interpolations
using DimensionalData: Dimension
using DimensionalData.Lookups
using DimensionalData.Dimensions: label
export shift, rebuild_with_slice, cat_arpes


"""
    shift(
        data::AbstractDimArray,
        shift_dim::Union{Symbol, Dimension},
        other_dim::Union{Symbol, Dimension},
        values::AbstractVector{<:Real},
    )

Shift the `data` array along `shift_dim` by an amount specified in `values`, for each position along `other_dim`.

- `data`: An `AbstractDimArray` to be shifted.
- `shift_dim`: The dimension along which to shift (as a `Symbol` or `Dimension`).
- `other_dim`: The dimension whose index determines which value from `values` to use.
- `values`: A vector of shift values, one for each index along `other_dim`.

Returns a new `DimArray` with the same shape and metadata as `data`, shifted accordingly.

Throws `ArgumentError` if the dimensions or value lengths do not match requirements.
"""
function shift(
    data::AbstractDimArray,
    shift_dim::Union{Symbol,Dimension},
    other_dim::Union{Symbol,Dimension},
    values::AbstractVector{<:Real},
)

    if length(values) != size(data, other_dim)
        error_msg = "Length of values vector ($(length(values))) does not match the size of the other dimension ($(size(data, other_dim)))."
        throw(ArgumentError(error_msg))
    end
    if !hasdim(data, shift_dim)
        error_msg = "Dimension :$shift_dim not found in the data object."
        throw(ArgumentError(error_msg))
    end
    if !hasdim(data, other_dim)
        error_msg = "Dimension :$other_dim not found in the data object."
        throw(ArgumentError(error_msg))
    end
    if ndims(data) < 2
        error_msg = "The dimension of the ARPESData object must be larger than two-dimensional for this operation."
        throw(ArgumentError(error_msg))
    end

    sdim = dimnum(data, shift_dim)
    odim = dimnum(data, other_dim)

    coords_vec = collect(dims(data, shift_dim))
    if !_is_equal_spacing(coords_vec)
        return _shift_irregular(data, shift_dim, values)
    end
    c1 = coords_vec[1]
    Δ = coords_vec[2] - c1
    itp = extrapolate(interpolate(data, BSpline(Linear()), OnGrid()), NaN)

    result = similar(data)
    inds = CartesianIndices(data)

    result .= (i -> begin
        # shift amount corresponding to this position along shift_dim

        # compute interpolated index
        target_idx = (coords_vec[i[sdim]] - values[i[odim]] - c1) / Δ + 1

        idx_tuple = Tuple(i)
        return itp(Base.setindex(idx_tuple, target_idx, sdim)...)
    end).(inds)

    return DimArray(result; dims = dims(data))
end

"""
    shift_optimized(data::AbstractDimArray,
                    shift_dim::Union{Symbol, Dimension},
                    value::Real)

Shift `data` along `shift_dim` by a scalar `value`.

- `shift_dim` must be equally spaced.
- A single shift value is applied uniformly along `shift_dim`.
- Subpixel shift is performed via linear interpolation.
- Returns an array with the same shape and metadata.

# Arguments
- `data` : AbstractDimArray
- `shift_dim` : Dimension to shift along
- `value` : scalar shift amount

# Returns
- AbstractDimArray with the same size and metadata
"""
function shift(data::AbstractDimArray, shift_dim::Union{Symbol,Dimension}, value::Real)

    sdim = dimnum(data, shift_dim)

    coords_vec = collect(dims(data, shift_dim))
    if !_is_equal_spacing(coords_vec)
        return _shift_irregular(data, shift_dim, value)
    end
    c1 = coords_vec[1]
    Δ = coords_vec[2] - c1

    itp = extrapolate(interpolate(data, BSpline(Linear()), OnGrid()), NaN)

    result = similar(data)
    inds = CartesianIndices(data)

    result .= (i -> begin
        # compute interpolated index with uniform shift
        target_idx = (coords_vec[i[sdim]] - value - c1) / Δ + 1
        idx_tuple = Tuple(i)
        itp(Base.setindex(idx_tuple, target_idx, sdim)...)
    end).(inds)

    return DimArray(result; dims = dims(data))
end

"""
    shift(
        data::AbstractDimArray,
        shift_dim::Union{Symbol, Dimension},
        values::AbstractDimArray,
    )

Shift an `AbstractDimArray` object along `shift_dim` by values provided in a one-dimensional `AbstractDimArray`.

- `data`: An `AbstractDimArray` object to be shifted.
- `shift_dim`: The dimension along which to shift (as a `Symbol` or `Dimension`).
- `values`: A one-dimensional `AbstractDimArray` of shift values.

Returns a new `DimArray` with the same shape and metadata as `data`, shifted accordingly.

Throws `ArgumentError` if `values` is not one-dimensional.
"""
function shift(
    data::AbstractDimArray,
    shift_dim::Union{Symbol,Dimension},
    values::AbstractDimArray,
)
    if ndims(values) != 1
        error_msg = "Values must be a one-dimensional array."
        throw(ArgumentError(error_msg))
    end
    other_dim = name(dims(data))[1]
    values = parent(values)
    shift(data, shift_dim, other_dim, values)
end

function _shift_irregular(
    data::AbstractDimArray,
    shift_dim::Union{Symbol,Dimension},
    values::AbstractVector{<:Real},
)

    sdim = dimnum(data, shift_dim)

    coords_vec = collect(dims(data, shift_dim))
    N = length(coords_vec)
    @assert length(values) == N

    # Check monotonicity (either strictly increasing or decreasing)
    if !(all(diff(coords_vec) .> 0) || all(diff(coords_vec) .< 0))
        throw(ArgumentError("shift_dim must be monotonic"))
    end

    # Build inverse mapping: coordinate → index via linear interpolation
    idxs = collect(1:N)
    coord_to_index = extrapolate(interpolate((coords_vec,), idxs, Gridded(Linear())), NaN)

    # Interpolator for the data itself (linear, on-grid)
    itp = extrapolate(interpolate(data, BSpline(Linear()), OnGrid()), NaN)

    result = similar(data)
    inds = CartesianIndices(data)

    result .= (i -> begin
        # Original physical coordinate along shift_dim
        x = coords_vec[i[sdim]]

        # Apply shift in physical units
        x_shifted = x - values[i[sdim]]

        # Convert shifted coordinate to fractional index (irregular grid)
        target_idx = coord_to_index(x_shifted)

        idx_tuple = Tuple(i)
        itp(Base.setindex(idx_tuple, target_idx, sdim)...)
    end).(inds)

    return DimArray(result; dims = dims(data))
end


"""
    _shift_irregular(data::AbstractDimArray,
                    shift_dim::Union{Symbol, Dimension},
                    value::Real)

Shift `data` along `shift_dim` by a scalar `value`, supporting irregular grids.

- `shift_dim` can be non-uniform but must be monotonic.
- A single shift value is applied uniformly.
- Subpixel shift is performed via linear interpolation.
- Returns an array with the same shape and metadata.

# Arguments
- `data` : AbstractDimArray
- `shift_dim` : Dimension to shift along
- `value` : scalar shift amount

# Returns
- AbstractDimArray with the same size and metadata
"""
function _shift_irregular(
    data::AbstractDimArray,
    shift_dim::Union{Symbol,Dimension},
    value::Real,
)

    sdim = dimnum(data, shift_dim)

    coords_vec = collect(dims(data, shift_dim))
    N = length(coords_vec)

    # Check monotonicity (required for inverse mapping)
    if !(all(diff(coords_vec) .> 0) || all(diff(coords_vec) .< 0))
        throw(ArgumentError("shift_dim must be monotonic"))
    end

    # Build coordinate → index mapping (inverse map)
    idxs = collect(1:N)
    coord_to_index = extrapolate(interpolate((coords_vec,), idxs, Gridded(Linear())), NaN)

    # Interpolator for data (index space)
    itp = extrapolate(interpolate(data, BSpline(Linear()), OnGrid()), NaN)

    result = similar(data)
    inds = CartesianIndices(data)

    result .= (i -> begin
        # Original coordinate
        x = coords_vec[i[sdim]]

        # Apply uniform shift
        x_shifted = x - value

        # Convert coordinate → fractional index
        target_idx = coord_to_index(x_shifted)

        idx_tuple = Tuple(i)
        itp(Base.setindex(idx_tuple, target_idx, sdim)...)
    end).(inds)

    return DimArray(result; dims = dims(data))
end

# --------------------------------------------

# Ensure `dimsel` fixes exactly one dimension (so A[dimsel] is (n-1)-D).
function _target_slice(A::AbstractDimArray, dimsel)
    fullN = ndims(A)
    target = A[dimsel]
    sliceN = ndims(target)

    sliceN == fullN - 1 || throw(
        ArgumentError(
            "dimsel must fix exactly one dimension: ndims(A)=$fullN, ndims(A[dimsel])=$sliceN",
        ),
    )

    return target
end

"""
    rebuild_with_slice(A::DimArray, dimsel, X::AbstractArray) -> DimArray

Return a new `DimArray` by copying `A` and replacing the (n-1)-D slice
specified by `dimsel` (e.g. `Ti(At(t0))`) with `X`.

Checks performed (always):
- `dimsel` must reduce dimensionality by exactly 1
- `size(A[dimsel]) == size(X)`
"""
function rebuild_with_slice(A::AbstractDimArray, dimsel, X::AbstractArray)
    target = _target_slice(A, dimsel)

    size(target) == size(X) ||
        throw(DimensionMismatch("size(A[dimsel])=$(size(target)) != size(X)=$(size(X))"))

    B = copy(A)
    B[dimsel] = X
    return B
end

"""
    rebuild_with_slice(A::DimArray, dimsel, X::DimArray) -> DimArray

Same as the `AbstractArray` method, but also checks dimension metadata.

Checks performed (always):
- `dimsel` must reduce dimensionality by exactly 1
- `size(A[dimsel]) == size(X)`
- `dims(A[dimsel]) == dims(X)`

Write policy:
- Writes `parent(X)` so that the returned array keeps `A`'s dims/metadata.
"""
function rebuild_with_slice(A::AbstractDimArray, dimsel, X::DimArray)
    target = _target_slice(A, dimsel)

    size(target) == size(X) ||
        throw(DimensionMismatch("size(A[dimsel])=$(size(target)) != size(X)=$(size(X))"))

    dims(target) == dims(X) ||
        throw(ArgumentError("dims(A[dimsel])=$(dims(target)) != dims(X)=$(dims(X))"))

    B = copy(A)
    B[dimsel] = parent(X)
    return B
end

# --------------------------------------------

"""
    cat_arpes(args::ARPESData...; dims::Union{Symbol,Dim})

Concatenate multiple `ARPESData` objects along the specified dimension.

# Arguments
- `args::ARPESData...`: One or more `ARPESData` objects to concatenate.
- `dims::Union{Symbol,Dim}`: The dimension (as a `Symbol` or `Dim`) along which to concatenate.

# Throws
- `ArgumentError` if any of the input data objects do not have the specified dimension.

# Returns
- A new `ARPESData` object resulting from concatenation along the specified dimension.

# Example
```julia
cat_arpes(data1, data2; dims=:energy)
```
"""
function cat_arpes(args::ARPESData...; dims::Union{Symbol,Dimension})
    missing_indices = findall(d -> !hasdim(d, dims), args)
    if !isempty(missing_indices)
        error_msg = "Dimension :$dims not found in the following argument(s): $(join(missing_indices, ", "))"
        throw(ArgumentError(error_msg))
    end
    lookups = map(d -> lookup(d, dims), args)
    combined_lookup_data = vcat(map(collect, lookups)...)
    proto_lookup = first(lookups)
    target_dims = Dim{dims}(combined_lookup_data, metadata = metadata(proto_lookup))
    return cat(args...; dims = target_dims)
end

_default_rtol(T) = T <: AbstractFloat ? √eps(T) : 0

function _is_equal_spacing(x::AbstractVector; atol = 0.0, rtol = _default_rtol(eltype(x)))
    if length(x) < 2
        return true
    end
    diffs = diff(x)
    return all(isapprox.(diff(diffs), 0; atol = atol, rtol = rtol))
end
