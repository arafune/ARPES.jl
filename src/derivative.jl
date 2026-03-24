using LinearAlgebra
using DimensionalData
using DimensionalData: rebuild
using DimensionalData: DimArray
using ..ARPES: _apply_along_dim

export derivative

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
  interior and forward/backward differences are used at the boundaries.
- For nonuniform coordinates, the first derivative is computed using the local
  spacing of neighboring points.
- For `order > 1` on nonuniform grids, higher-order derivatives are obtained by
  repeated application of the first-derivative operator, which is an
  approximation.
- When `dim` does not exist in `A`, the corresponding error from
  `DimensionalData.dims` is propagated.
- `order < 0` throws an error.
"""
function derivative(A::AbstractDimArray, dim; order::Integer = 1)
    # NOTE: higher-order derivatives for nonuniform grids
    # are computed by repeated first derivatives (approximate)
    order < 0 && error("order must be ≥ 0")
    order == 0 && return A

    d = DimensionalData.dims(A, dim)
    x = collect(d)

    if _is_equal_spacing(x)
        Δx = x[2] - x[1]
        f = y -> _diff_uniform(y, Δx, order)
    else
        f = y -> _diff_nonuniform(y, x, order)
    end

    return _apply_along_dim(A, dim, f)
end

"""
    _diff_uniform(y::AbstractVector, Δx::Real, order::Int)

Compute a numerical derivative of order `order` for uniformly spaced samples
`y` with spacing `Δx`.

This helper promotes the input to floating-point values and applies the
first-order finite-difference operator repeatedly.
"""
function _diff_uniform(y::AbstractVector, Δx::Real, order::Int)
    out = copy(float.(y))

    for _ = 1:order
        out = _central_diff_uniform(out, Δx)
    end

    return out
end

"""
    _central_diff_uniform(y, Δx)

Compute the first derivative of uniformly spaced samples `y` using central
differences in the interior and first-order one-sided differences at the
boundaries.
"""
function _central_diff_uniform(y, Δx)
    n = length(y)
    out = similar(y)

    # central difference (interior)
    @inbounds for i = 2:(n-1)
        out[i] = (y[i+1] - y[i-1]) / (2Δx)
    end

    # forward/backward at boundaries
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
array with the original metadata preserved.
"""
function minimum_gradient(
    A::AbstractDimArray,
    dims::Tuple{T,T},
) where {T<:Union{DimensionalData.Dimension,Symbol}}
    gradients = A ./ _gradient_modulus(A)
    return rebuild(A, data = gradients)
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
    maximum_curvature(A::AbstractDimArray, dim, alpha::Real=0.1)

Compute the 1D maximum-curvature transform of `A` along `dim`.

This method uses the first and second derivatives along `dim` to enhance
dispersive features. The tuning parameter `alpha` must be positive and controls
the regularization strength in the curvature denominator.
"""
function maximum_curvature(
    A::AbstractDimArray,
    dim::T,
    alpha::Real = 0.1,
) where {T<:Union{DimensionalData.Dimension,Symbol}}
    @assert alpha > 0
    dim = DimensionalData.dims(A, dim)
    dx = _step(dim)
    df = derivative(A, dim)
    d2f = derivative(A, dim, order = 2)
    denominator = (alpha * abs(minimum(df)) ^ 2 .+ d2f .^ 2) .^ 1.5
    curv = d2f ./ denominator
    return rebuild(A, data = curv)
end

"""
    maximum_curvature(A::AbstractDimArray, dims::Tuple, alpha::Real=0.1, weight2d::Real=1.0)

Compute the 2D maximum-curvature transform of `A` over the dimension pair
`dims`.

First- and second-order derivatives, including the mixed derivative, are used
to build the 2D curvature expression. `alpha` must be positive, and `weight2d`
controls the relative scaling between the two dimensions.
"""
function maximum_curvature(
    A::AbstractDimArray,
    dims::Tuple{T,T},
    alpha::Real = 0.1,
    weight2d::Real = 1.0,
) where {T<:Union{DimensionalData.Dimension,Symbol}}
    @assert weight2d != 0
    @assert alpha > 0
    dim1 = DimensionalData.dims(A, dims[1])
    dim2 = DimensionalData.dims(A, dims[2])
    dx, dy = _step(dim1), _step(dim2)
    df = derivative(A, dim1), derivative(A, dim2)
    d2f = (
        derivative(A, dim1, order = 2),
        derivative(A, dim2, order = 2),
        derivative(df[1], dim2),
    )
    weight = weight2d > 0 ? (dx / dy)^2 * weight2d : (dx / dy)^2 / abs(weight2d)
    avg_x = abs(minimum(df[1]))
    avg_y = abs(minimum(df[2]))
    avg = maximum(avg_x^2, weight * avg_y^2)
    numerator =
        (alpha * avg + weight * df[1] .* df[2]) .* d2f[2] -
        2 * weight * df[1] .* df[2] .* d2f[3] +
        weight * (alpha * avg + df[2] .* df[2]) .* d2f[1]
    denominator = (alpha * avg + weight * df[1] .^ 2 + df[2] .^ 2) .^ 1.5
    curv = numerator ./ denominator
    return rebuild(A, data = curv)
end
