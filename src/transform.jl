using Interpolations
using DimensionalData: Dimension, Bins
using DimensionalData.Lookups
using DimensionalData.Dimensions: label
export rebin, shift, rebuild_with_slice



"""
    rebin(data::AbstractDimArray, dim::Union{Symbol,Dimension}, bins::Int)

Rebins the data along the specified dimension `dim` into `bins` number of bins.
This is a convenience wrapper that constructs a `Bins` object from the integer
and calls the main `rebin` method.

# Arguments
- `data`: The input `AbstractDimArray` to be rebinned.
- `dim`: The dimension (as a `Symbol` or `Dimension`) along which to rebin.
- `bins`: The number of bins to use for rebinning.

# Returns
A new array with the specified dimension rebinned into the given number of bins.
"""
function rebin(data::AbstractDimArray, dim::Union{Symbol,Dimension}, bins::Int)
    bins = Bins(bins)
    rebin(data, dim, bins)
end

"""
    rebin(data::AbstractDimArray, dim::Union{Symbol,Dimension}, bins::Bins)

Rebins the data along the specified dimension `dim` using the provided `Bins` object.
Groups the data along `dim` into bins, computes the mean within each bin, and concatenates
the results along the rebinned dimension.

# Arguments
- `data`: The input `AbstractDimArray` to be rebinned.
- `dim`: The dimension (as a `Symbol` or `Dimension`) along which to rebin.
- `bins`: A `Bins` object specifying the binning strategy.

# Returns
A new array with the specified dimension rebinned according to `bins`.
"""
function rebin(data::AbstractDimArray, dim::Union{Symbol,Dimension}, bins::Bins)
    @assert hasdim(data, dim) "Dimension :$dim not found in the data object."
    dim = dims(data, dim)
    groupby_dim_array = groupby(data, dim=>bins)
    binned_data = map(a -> mean(a, dims = dim), groupby_dim_array)
    reduce((x, y) -> cat(x, y; dims = name(dim)), binned_data)
end

"""
    shift(
        data::AbstractDimArray,
        shift_dim::Union{Symbol, Dimension},
        other_dim::Union{Symbol, Dimension},
        values::AbstractVector{<:Real},
    )

Shift the `data` array along `shift_dim` by an amount specified in `values`,
for each position along `other_dim`.

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

    if ndims(data) < 2
        error_msg = "The dimension of the $data object must be larger than two-dimensional for this operation."
        throw(ArgumentError(error_msg))
    end
    if !hasdim(data, other_dim)
        error_msg = "Dimension :$other_dim not found in the data object."
        throw(ArgumentError(error_msg))
    end
    if length(values) != size(data, other_dim)
        error_msg = "Length of values vector ($(length(values))) does not match the size of the other dimension ($(size(data, other_dim)))."
        throw(ArgumentError(error_msg))
    end
    if !hasdim(data, shift_dim)
        error_msg = "Dimension :$shift_dim not found in the data object."
        throw(ArgumentError(error_msg))
    end

    sdim = dimnum(data, shift_dim)
    odim = dimnum(data, other_dim)

    coords_vec = collect(dims(data, shift_dim))
    if !_is_equal_spacing(coords_vec)
        return _shift_irregular(data, shift_dim, other_dim, values)
    end
    c1 = coords_vec[1]
    Δ = coords_vec[2] - c1
    itp = extrapolate(interpolate(data, BSpline(Linear())), NaN)

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

    itp = extrapolate(interpolate(data, BSpline(Linear())), NaN)

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
    @debug "values are AbstractDimArray"
    if ndims(values) != 1
        error_msg = "Values must be a one-dimensional array."
        throw(ArgumentError(error_msg))
    end
    other_dim = name(dims(values))[1]
    values = parent(values)
    shift(data, shift_dim, other_dim, values)
end

function _shift_irregular(
    data::AbstractDimArray,
    shift_dim::Union{Symbol,Dimension},
    other_dim::Union{Symbol,Dimension},
    values::AbstractVector{<:Real},
)
    sdim = dimnum(data, shift_dim)
    odim = dimnum(data, other_dim)

    coords = collect(parent(dims(data, shift_dim)))
    if !(all(diff(coords) .> 0) || all(diff(coords) .< 0))
        throw(ArgumentError("shift_dim must be monotonic"))
    end
    A = parent(data)

    result = similar(A)

    inds = CartesianIndices(A)

    result .= (i -> begin
        idx = Tuple(i)

        shift_val = values[i[odim]]

        slice = view(A, ntuple(d -> d == sdim ? Colon() : idx[d], ndims(A))...)

        itp = extrapolate(interpolate((coords,), slice, Gridded(Linear())), NaN)

        x = coords[i[sdim]]
        x_shifted = x - shift_val

        return itp(x_shifted)
    end).(
        inds,
    )

    return DimArray(result, dims(data))
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
    itp = extrapolate(interpolate(data, BSpline(Linear())), NaN)

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
function _get_validated_slice(A, dimsel)
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
    rebuild_with_slice(A::AbstractDimArray, dimsel, X::AbstractArray) -> DimArray

Return a new `DimArray` by copying `A` and replacing the (n-1)-D slice
specified by `dimsel` (e.g. `Ti(At(t0))`) with `X`.

# NOTE:
# Only values of `X` are written. The resulting array keeps `A`'s dims/metadata.

- `dimsel` must fix exactly one dimension (e.g. `Ti(At(t0))`),
  leaving all other dimensions unmodified.

Checks performed (always):
- `dimsel` must fix exactly one dimension
- `size(A[dimsel]) == size(X)`
"""
function rebuild_with_slice(A::AbstractDimArray, dimsel, X::AbstractArray)
    target = _get_validated_slice(A, dimsel)
    size(target) == size(X) ||
        throw(DimensionMismatch("size(A[dimsel])=$(size(target)) != size(X)=$(size(X))"))

    B = copy(A)
    B[dimsel] = X
    return B
end

"""
    rebuild_with_slice(A::AbstractDimArray, dimsel, X::AbstractDimArray) -> DimArray

Same as the `AbstractArray` method, but also checks dimension metadata.

Checks performed (always):
- `size(A[dimsel]) == size(X)`
- `dims(A[dimsel]) == dims(X)`
- `dimsel` must fix exactly one dimension (e.g. `Ti(At(t0))`)
- Selecting multiple indices or no indices is not allowed

Write policy:
- Writes `parent(X)` so that the returned array keeps `A`'s dims/metadata.
"""
function rebuild_with_slice(A::AbstractDimArray, dimsel, X::AbstractDimArray)
    target = _get_validated_slice(A, dimsel)

    size(target) == size(X) ||
        throw(DimensionMismatch("size(A[dimsel])=$(size(target)) != size(X)=$(size(X))"))

    dims(target) == dims(X) ||
        throw(ArgumentError("dims(A[dimsel])=$(dims(target)) != dims(X)=$(dims(X))"))

    B = copy(A)
    B[dimsel] = parent(X)
    return B
end


