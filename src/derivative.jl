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
- `A::AbstractDimArray`: Input data.
- `dim`: Dimension along which the derivative is computed. This is resolved with
  `DimensionalData.dims(A, dim)`, so it can be given in the same forms accepted
  there (for example a dimension label or dimension object).
- `order::Integer=1`: Derivative order. `order = 0` returns `A` unchanged.

# Returns
- A new `AbstractDimArray` with the same dimensions as `A`, containing the
  numerical derivative along `dim`.

# Notes
- For equally spaced coordinates, a central-difference stencil is used in the
  interior and forward/backward differences are used at the boundaries.
- For nonuniform coordinates, the first derivative is computed using the local
  spacing of neighboring points.
- For `order > 1` on nonuniform grids, higher-order derivatives are obtained by
  repeated application of the first-derivative operator, which is an
  approximation.
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

function _diff_uniform(y::AbstractVector, Δx::Real, order::Int)
    out = copy(float.(y))

    for _ = 1:order
        out = _central_diff_uniform(out, Δx)
    end

    return out
end

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

function _diff_nonuniform(y::AbstractVector, x::AbstractVector, order::Int)
    out = copy(float.(y))

    for _ = 1:order
        out = _central_diff_nonuniform(out, x)
    end

    return out
end

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


function minimum_gradient(
    A::AbstractDimArray,
    dims::Tuple{T,T},
) where {T<:Union{DimensionalData.Dimension,Symbol}}
    gradients = A ./ _gradient_modulus(A)
    return rebuild(A, data = gradients)
end

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
