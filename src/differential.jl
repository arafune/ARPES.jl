using DimensionalData
using ..ARPES: _apply_along_dim


function differential(A::AbstractDimArray, dim, order::Integer)
    order < 0 && error("order must be ≥ 0")
    order == 0 && return A

    d = DimensionalData.dims(A, dim)
    x = collect(d)

    if _is_equal_spacing(x)
        Δx = float(x[2] - x[1])
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
