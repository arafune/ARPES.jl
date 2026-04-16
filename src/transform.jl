using Interpolations
using DimensionalData: Bins
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
    shift(data::AbstractDimArray, shift_dim::Union{Symbol, Dimension}, other_dim::Union{Symbol, Dimension}, values::AbstractVector{<:Real})
    shift(data::AbstractDimArray, shift_dim::Union{Symbol, Dimension}, values::AbstractDimArray)
    shift(data::AbstractDimArray, shift_dim::Union{Symbol, Dimension}, value::Real)

Shift `data` along `shift_dim` using linear interpolation.

# Arguments
- `data`: The input `AbstractDimArray` to be shifted.
- `shift_dim`: The dimension (as a `Symbol` or `Dimension`) along which to shift.

- `value`: Apply a uniform scalar shift along `shift_dim`.
- `other_dim`, `values`: Apply shifts that vary along `other_dim`, using one value for each index along that dimension.
- `values::AbstractDimArray`: A one-dimensional shorthand where the dimension of `values` determines `other_dim`.
- `dim_extend=false`: When `true`, rebuild `shift_dim` on an extended axis with the same step size.

For equally spaced coordinates, shifting is performed directly in index space.
For monotonic irregular coordinates, the implementation falls back to irregular-grid handling.
Dimension extension is supported only for equally spaced coordinates.

# Returns
A new `DimArray` with the shifted data. When `dim_extend=true`, the returned array is larger along `shift_dim`.

# Throws
- `ArgumentError`: If `values` is not one-dimensional, if a requested dimension is missing, if the length of `values` does not match the size of `other_dim`, or if `dim_extend=true` is used with a non-equally-spaced `shift_dim`.
"""
function shift(
    data::AbstractDimArray,
    shift_dim::Union{Symbol,Dimension},
    other_dim::Union{Symbol,Dimension},
    values::AbstractVector{<:Real},
    ;
    dim_extend::Bool = false,
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
    if dim_extend
        _validate_uniform_shift_coords(coords_vec)
        return _shift_with_dim_extend(data, shift_dim, other_dim, values)
    end
    if !_is_equal_spacing(coords_vec)
        return _shift_irregular(data, shift_dim, other_dim, values)
    end
    c1, Δ = _validate_uniform_shift_coords(coords_vec)
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

function shift(
    data::AbstractDimArray,
    shift_dim::Union{Symbol,Dimension},
    value::Real;
    dim_extend::Bool = false,
)

    sdim = dimnum(data, shift_dim)

    coords_vec = collect(dims(data, shift_dim))
    if dim_extend
        _validate_uniform_shift_coords(coords_vec)
        return _shift_with_dim_extend(data, shift_dim, value)
    end
    if !_is_equal_spacing(coords_vec)
        return _shift_irregular(data, shift_dim, value)
    end
    c1, Δ = _validate_uniform_shift_coords(coords_vec)

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

function shift(
    data::AbstractDimArray,
    shift_dim::Union{Symbol,Dimension},
    values::AbstractDimArray,
    ;
    dim_extend::Bool = false,
)
    @debug "values are AbstractDimArray"
    other_dim = only(name(dims(values)))
    values = parent(values)
    shift(data, shift_dim, other_dim, values; dim_extend)
end

function _build_extended_shift_dim(dim::Dimension, shift_range::Tuple{<:Real,<:Real})
    coords = collect(parent(lookup(dim)))
    c1, Δ = _validate_uniform_shift_coords(coords)
    min_shift, max_shift = shift_range
    extension_steps = ceil(Int, (max_shift - min_shift) / abs(Δ))
    start = c1 + (Δ > 0 ? min_shift : max_shift)
    extended_lookup = range(start, step = Δ, length = length(coords) + extension_steps)
    return rebuild(dim; val = extended_lookup)
end

function _shift_with_dim_extend(
    data::AbstractDimArray,
    shift_dim::Union{Symbol,Dimension},
    other_dim::Union{Symbol,Dimension},
    values::AbstractVector{<:Real},
)
    sdim = dimnum(data, shift_dim)
    odim = dimnum(data, other_dim)
    shift_axis = dims(data, shift_dim)
    coords_vec = collect(parent(lookup(shift_axis)))
    c1, Δ = _validate_uniform_shift_coords(coords_vec)

    extended_dim = _build_extended_shift_dim(shift_axis, (minimum(values), maximum(values)))
    extended_coords = collect(parent(lookup(extended_dim)))
    out_size = ntuple(d -> d == sdim ? length(extended_coords) : size(data, d), ndims(data))
    result = similar(parent(data), out_size)
    itp = extrapolate(interpolate(data, BSpline(Linear())), NaN)

    result .= (i -> begin
        target_idx = (extended_coords[i[sdim]] - values[i[odim]] - c1) / Δ + 1
        idx_tuple = Tuple(i)
        return itp(Base.setindex(idx_tuple, target_idx, sdim)...)
    end).(
        CartesianIndices(result),
    )

    new_dims = Base.setindex(dims(data), extended_dim, sdim)
    return rebuild(data; data = result, dims = new_dims)
end

function _shift_with_dim_extend(
    data::AbstractDimArray,
    shift_dim::Union{Symbol,Dimension},
    value::Real,
)
    sdim = dimnum(data, shift_dim)
    shift_axis = dims(data, shift_dim)
    coords_vec = collect(parent(lookup(shift_axis)))
    c1, Δ = _validate_uniform_shift_coords(coords_vec)

    extended_dim = _build_extended_shift_dim(shift_axis, (value, value))
    extended_coords = collect(parent(lookup(extended_dim)))
    out_size = ntuple(d -> d == sdim ? length(extended_coords) : size(data, d), ndims(data))
    result = similar(parent(data), out_size)
    itp = extrapolate(interpolate(data, BSpline(Linear())), NaN)

    result .= (i -> begin
        target_idx = (extended_coords[i[sdim]] - value - c1) / Δ + 1
        idx_tuple = Tuple(i)
        return itp(Base.setindex(idx_tuple, target_idx, sdim)...)
    end).(
        CartesianIndices(result),
    )

    new_dims = Base.setindex(dims(data), extended_dim, sdim)
    return rebuild(data; data = result, dims = new_dims)
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

function _validate_uniform_shift_coords(coords_vec::AbstractVector)
    length(coords_vec) >= 2 ||
        throw(ArgumentError("shift_dim must contain at least two coordinate values."))
    _is_equal_spacing(coords_vec) ||
        throw(ArgumentError("shift_dim must be equally spaced for this operation."))
    Δ = coords_vec[2] - coords_vec[1]
    iszero(Δ) &&
        throw(ArgumentError("shift_dim must have non-zero spacing for this operation."))
    return coords_vec[1], Δ
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
