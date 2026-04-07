using LinearAlgebra
using DimensionalData
using DimensionalData: rebuild
using DimensionalData: DimArray
using ..ARPES: _apply_along_dim, _apply_along_dims

export derivative, curvature, minimum_gradient

"""
    derivative(A::AbstractDimArray, dim; order::Integer=1)

Compute the numerical derivative of `A` along the dimension `dim`.

The derivative is evaluated independently along each 1D line selected by `dim`,
while preserving the original array shape, dimensions, and metadata.

# Arguments
- `A::AbstractDimArray`: Input data. Each slice along `dim` is differentiated
  independently.
- `dim`: Dimension along which the derivative is computed. This is resolved with
  `DimensionalData.dims(A, dim)`, so it can be given in the same forms accepted
  there (for example a dimension label or dimension object).
- `order::Integer=1`: Derivative order. `order = 0` returns `A` unchanged.

# Returns
- A new `AbstractDimArray` with the same dimensions as `A`, containing the
  numerical derivative along `dim`. The computed values are promoted to a
  floating-point representation.

# Requirements
- For `order > 0`, the lookup of `dim` must contain at least two coordinate
  values.
- When the coordinates are treated as equally spaced, the current implementation
  uses second-order one-sided boundary stencils and therefore requires at least
  three samples along `dim`.

# Examples
```julia
julia> using ARPES, DimensionalData

julia> t = range(0, 2π; length = 200);

julia> A = DimArray(sin.(t), (Dim{:t}(t),));

julia> dA = derivative(A, :t);

julia> maximum(abs.(parent(dA) .- cos.(t))) < 1e-2
true
```

Different slices of a higher-dimensional array are processed independently:

```julia
julia> data = rand(50, 10);

julia> A2 = DimArray(data, (Dim{:eV}(range(-1, 1; length = 50)), Dim{:phi}(1:10)));

julia> dA2 = derivative(A2, :eV);

julia> size(dA2) == size(A2)
true
```

# Notes
- For equally spaced coordinates, a central-difference stencil is used in the
  interior and, by default, second-order one-sided differences are used at the
  boundaries.
- For nonuniform coordinates, the implementation currently uses a three-point
  formula in the interior and first-order forward/backward differences at the
  boundaries.
- For `order > 1` on nonuniform grids, higher-order derivatives are obtained by
  repeated application of that first-derivative operator. In other words, the
  nonuniform-grid path is currently based on repeated first-derivative
  approximations rather than a dedicated higher-order scheme.
- When `dim` does not exist in `A`, the corresponding error from
  `DimensionalData.dims` is propagated.
- `order < 0` throws an error.
"""
function derivative(A::AbstractDimArray, dim; order::Integer = 1, accuracy::Int = 2)
    # NOTE: higher-order derivatives for nonuniform grids
    # are computed by repeated first derivatives (approximate)
    order < 0 && error("order must be ≥ 0")
    order == 0 && return A

    d = DimensionalData.dims(A, dim)
    x = collect(d)

    if _is_equal_spacing(x)
        Δx = x[2] - x[1]
        f = y -> _diff_uniform(y, Δx, order; accuracy = accuracy)
    else
        f = y -> _diff_nonuniform(y, x, order)
    end

    return _apply_along_dim(A, dim, f)
end

"""
    _diff_uniform(y::AbstractVector, Δx::Real, order::Int; accuracy::Int=2)

Compute a numerical derivative of order `order` for uniformly spaced samples
`y` with spacing `Δx`.

This helper promotes the input to floating-point values and applies the
first-order finite-difference operator repeatedly. The keyword `accuracy`
selects the boundary stencil order passed to `_central_diff_uniform`.
"""
function _diff_uniform(y::AbstractVector, Δx::Real, order::Int; accuracy::Int = 2)
    out = copy(float.(y))

    for _ = 1:order
        out = _central_diff_uniform(out, Δx; accuracy = accuracy)
    end

    return out
end

"""
    _central_diff_uniform(y, Δx; accuracy::Int=2)

Compute the first derivative of uniformly spaced samples `y` using central
differences in the interior.

The keyword `accuracy` controls the one-sided stencil used at the boundaries:

- `accuracy = 2` (default): second-order one-sided differences
- `accuracy = 1`: first-order forward/backward differences
"""

_central_diff_uniform(y, Δx; accuracy::Int = 2) =
    _central_diff_uniform(y, Δx, Val(accuracy))
function _central_diff_uniform(y, Δx, ::Val{2})
    n = length(y)
    out = similar(y)
    @inbounds for i = 2:(n-1)
        out[i] = (y[i+1] - y[i-1]) / (2Δx)
    end
    out[1] = (-3y[1] + 4y[2] - y[3]) / (2Δx)
    out[end] = (3y[end] - 4y[end-1] + y[end-2]) / (2Δx)
    return out
end

function _central_diff_uniform(y, Δx, ::Val{1})
    n = length(y)
    out = similar(y)
    @inbounds for i = 2:(n-1)
        out[i] = (y[i+1] - y[i-1]) / (2Δx)
    end
    out[1] = (y[2] - y[1]) / Δx
    out[end] = (y[end] - y[end-1]) / Δx
    return out
end

"""
    _diff_nonuniform(y::AbstractVector, x::AbstractVector, order::Int)

Compute a numerical derivative of order `order` for samples `y` observed at the
nonuniform coordinates `x`.

Higher-order derivatives are constructed by repeated application of the
first-derivative operator implemented in `_central_diff_nonuniform`.
"""
function _diff_nonuniform(y::AbstractVector, x::AbstractVector, order::Int)
    out = copy(float.(y))

    for _ = 1:order
        out = _central_diff_nonuniform(out, x)
    end

    return out
end

"""
    _central_diff_nonuniform(y, x)

Compute the first derivative of samples `y` on the nonuniform coordinate grid
`x`.

Interior points use the three-point nonuniform finite-difference formula, while
the boundaries use first-order forward and backward differences.
"""
function _central_diff_nonuniform(y, x)
    n = length(y)
    out = similar(y)

    @inbounds for i = 2:(n-1)
        dx1 = x[i] - x[i-1]
        dx2 = x[i+1] - x[i]

        out[i] = (
            -dx2/(dx1*(dx1+dx2)) * y[i-1] +
            (dx2-dx1)/(dx1*dx2) * y[i] +
            dx1/(dx2*(dx1+dx2)) * y[i+1]
        )
    end

    # forward/backward (1st order)
    out[1] = (y[2] - y[1]) / (x[2] - x[1])
    out[end] = (y[end] - y[end-1]) / (x[end] - x[end-1])

    return out
end


"""
    minimum_gradient(A::AbstractDimArray, dims::Tuple)

Normalize `A` by the local gradient modulus computed over the specified pair of
dimensions.

This operation returns `A ./ _gradient_modulus(A)` rebuilt as a dimensional
array with the original metadata preserved. For arrays with additional
dimensions, the normalization is applied independently to each block selected by
`dims`.
"""
function minimum_gradient(
    A::AbstractDimArray,
    dims::Tuple{T,T},
) where {T<:Union{DimensionalData.Dimension,Symbol}}
    return _apply_along_dims(A, dims, _minimum_gradient_2d)
end

function _minimum_gradient_2d(A::AbstractDimArray)
    return parent(A) ./ _gradient_modulus(A)
end

"""
    _vector_diff(arr::AbstractArray, delta::NTuple{N,Int}, n::Int=1) where {N}

Compute finite differences of `arr` between slices separated by `delta`.

The tuple `delta` specifies the directional offset for each dimension. When
`n > 1`, the same directional difference is applied recursively `n` times.
"""
function _vector_diff(arr::AbstractArray, delta::NTuple{N,Int}, n::Int = 1) where {N}
    n == 0 && return arr
    @assert n > 0

    inds1 = ntuple(i -> begin
        d = delta[i]
        if d > 0
            (1+d):size(arr, i)
        elseif d < 0
            1:(size(arr, i)+d)
        else
            Colon()
        end
    end, N)

    inds2 = ntuple(i -> begin
        d = delta[i]
        if d > 0
            1:(size(arr, i)-d)
        elseif d < 0
            (1-d):size(arr, i)
        else
            Colon()
        end
    end, N)

    diff = arr[inds1...] .- arr[inds2...]

    return n > 1 ? _vector_diff(diff, delta, n - 1) : diff
end

"""
    _gradient_modulus(data::Union{AbstractDimArray,AbstractArray}; delta::Int=1)

Estimate the 2D gradient magnitude of `data` from directional finite
differences.

Eight forward, backward, and diagonal difference directions are evaluated, and
their Euclidean norm is returned for each point.
"""
function _gradient_modulus(data::Union{AbstractDimArray,AbstractArray}; delta::Int = 1)
    values = parent(data)
    grad = zeros(eltype(values), (8, size(values)...))
    grad[1, 1:(end-delta), :] .= _vector_diff(values, (delta, 0))
    grad[2, :, 1:(end-delta)] .= _vector_diff(values, (0, delta))
    grad[3, (delta+1):end, :] .= _vector_diff(values, (-delta, 0))
    grad[4, 1:(end-delta), 1:(end-delta)] .= _vector_diff(values, (delta, delta))
    grad[5, :, (delta+1):end] .= _vector_diff(values, (0, -delta))
    grad[6, 1:(end-delta), (delta+1):end] .= _vector_diff(values, (delta, -delta))
    grad[7, (delta+1):end, 1:(end-delta)] .= _vector_diff(values, (-delta, delta))
    grad[8, (delta+1):end, (delta+1):end] .= _vector_diff(values, (-delta, -delta))
    return dropdims(mapslices(norm, grad, dims = 1), dims = 1)
end

"""
    curvature(A::AbstractDimArray, dim, alpha::Real=0.1)

Compute the 1D maximum-curvature transform of `A` along `dim`.

This method uses the first and second derivatives along `dim` to enhance
dispersive features. The tuning parameter `alpha` must be positive and controls
the regularization strength in the curvature denominator. For arrays with
additional dimensions, the transform is applied independently to each slice
along `dim`.
"""
function curvature(
    A::AbstractDimArray,
    dim::T,
    alpha::Real = 0.1,
) where {T<:Union{DimensionalData.Dimension,Symbol}}
    @assert alpha > 0
    dim = DimensionalData.dims(A, dim)
    df, d2f = derivative(A, dim), derivative(A, dim, order = 2)
    demoninator = (alpha * maximum(abs.(parent(df)))^2 .+ parent(df) .^ 2) .^ 1.5
    return rebuild(A, parent(d2f) ./ demoninator)
end

"""
    curvature(A::AbstractDimArray, dims::Tuple, alpha::Real=0.1, weight2d::Real=1.0)

Compute the 2D curvature transform of `A` over the dimension pair
`dims`.

First- and second-order derivatives, including the mixed derivative, are used
to build the 2D curvature expression. `alpha` must be positive, and `weight2d`
controls the relative scaling between the two dimensions. For arrays with
additional dimensions, the transform is applied independently to each block
selected by `dims`.
"""
function curvature(
    A::AbstractDimArray,
    dims::Tuple{T,T},
    alpha::Real = 0.1,
    weight2d::Real = 1.0,
) where {T<:Union{DimensionalData.Dimension,Symbol}}
    @assert alpha > 0
    @assert weight2d != 0

    dim1 = DimensionalData.dims(A, dims[1])
    dim2 = DimensionalData.dims(A, dims[2])
    df = derivative(A, dim1), derivative(A, dim2)
    df1 = parent(df[1])
    df2 = parent(df[2])

    dx, dy = _step(collect(dim1)), _step(collect(dim2))
    weight = weight2d > 0 ? (dx / dy)^2 * weight2d : (dx / dy)^2 / abs(weight2d)
    scale_x = maximum(abs.(df1))^2
    scale_y = maximum(abs.(df2))^2
    avg_global = max(scale_x, weight * scale_y)

    return _apply_along_dims(
        A,
        dims,
        block -> _curvature_2d(block, dims, alpha, weight2d, avg_global),
    )
end

function _curvature_2d(
    A::AbstractDimArray,
    dims::Tuple{T,T},
    alpha::Real,
    weight2d::Real,
    avg::Real,
) where {T<:Union{DimensionalData.Dimension,Symbol}}
    dim1 = DimensionalData.dims(A, dims[1])
    dim2 = DimensionalData.dims(A, dims[2])
    dx, dy = _step(collect(dim1)), _step(collect(dim2))
    df = derivative(A, dim1), derivative(A, dim2)
    d2f = (
        derivative(A, dim1, order = 2),
        derivative(A, dim2, order = 2),
        derivative(df[1], dim2),
    )
    df1 = parent(df[1])
    df2 = parent(df[2])
    d2f1 = parent(d2f[1])
    d2f2 = parent(d2f[2])
    d2f12 = parent(d2f[3])

    weight = weight2d > 0 ? (dx / dy)^2 * weight2d : (dx / dy)^2 / abs(weight2d)

    numerator =
        (alpha * avg .+ weight .* df1 .* df2) .* d2f2 - 2 .* weight .* df1 .* df2 .* d2f12 +
        weight .* (alpha * avg .+ df2 .* df2) .* d2f1
    denominator = (alpha * avg .+ weight .* df1 .^ 2 .+ df2 .^ 2) .^ 1.5
    return numerator ./ denominator
end


