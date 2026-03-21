using Statistics
using DimensionalData
using LowLevelParticleFilters
using LinearAlgebra

export kalman_smooth_dim_llpf,
    boxcar_filter_dim, sg_filter_dim, binomial_filter_dim, gaussian_filter_dim
"""
    kalman_smooth_dim_llpf(A::AbstractDimArray, dim;
                           R=1.0, q=1.0, model=:rw)

Apply a 1D Kalman filter (using LowLevelParticleFilters.jl) along a specified
dimension of an N-dimensional DimArray.

The filter is applied independently to each slice along the chosen dimension,
i.e., this performs a broadcast/map of a 1D Kalman filter over the remaining
dimensions (no coupling between them).

Features:
- Supports non-uniform sampling (time intervals Δt are explicitly used)
- Handles missing observations (NaN values are skipped in the update step)

# Arguments
- A::AbstractDimArray:
    Input N-dimensional data.
- dim:
    Dimension (Symbol or Dimension) along which the filter is applied.
- R::Real:
    Observation noise variance.
- q::Real:
    Process noise intensity (controls smoothness).
- model::Symbol:
    State model to use:
    - `:rw` : Random walk (strong smoothing)
    - `:cv` : Constant velocity (smoother, allows linear trends)

# Returns
- DimArray:
    A new DimArray with the same shape as `A`, containing the filtered result.

# Notes
- This is not a fully coupled N-dimensional Kalman filter.
  Each slice along `dim` is processed independently.
- For a true N-dimensional Kalman filter, a joint state-space model
  (including spatial correlations) must be defined.
"""
function kalman_smooth_dim_llpf(A::AbstractDimArray, dim; R = 1.0, q = 1.0, model = :rw)

    d = DimensionalData.dims(A, dim)
    ax = collect(d)
    n = length(ax)

    # Δt for non-uniform sampling
    Δt = zeros(Float64, n)
    for i = 2:n
        Δt[i] = float(ax[i] - ax[i-1])
    end

    # State dimension
    nx = model == :rw ? 1 : model == :cv ? 2 : error("model must be :rw or :cv")

    out = similar(parent(A), eltype(A))
    dim_index = DimensionalData.dimnum(A, d)

    # Iterate over all slices except the target dimension
    other_sizes = ntuple(i -> i == dim_index ? 1 : size(A, i), ndims(A))

    for I in CartesianIndices(other_sizes)
        idx = ntuple(i -> i == dim_index ? Colon() : I[i], ndims(A))
        y = Array(A[idx...])

        # Initial state
        x̂ = zeros(nx)
        P = Matrix(I, nx, nx) * 1e3

        x_est = zeros(nx, n)

        for t = 1:n
            dt = Δt[t]

            if model == :rw
                F = [1.0]
                Q = [q * max(dt, 1e-12)]
                H = [1.0]
            else # :cv
                F = [
                    1.0 dt;
                    0.0 1.0
                ]
                Q = q * [
                    dt^3/3 dt^2/2;
                    dt^2/2 dt
                ]
                H = [1.0 0.0]
            end

            # Predict
            x̂ = F * x̂
            P = F * P * F' + Q

            # Update (skip NaN)
            if !isnan(y[t])
                S = H * P * H' + R
                K = P * H' / S
                x̂ = x̂ + K * (y[t] - H * x̂)
                P = (I - K * H) * P
            end

            x_est[:, t] = x̂
        end

        # Write back (use first state component)
        out[idx...] = x_est[1, :]
    end

    return rebuild(A, out)
end


# -------------------------------
# Utility: apply 1D function along a dimension
# -------------------------------
function _apply_along_dim(A::AbstractDimArray, dim, f::Function)
    d = DimensionalData.dims(A, dim)
    dim_index = DimensionalData.dimnum(A, d)

    out = similar(parent(A), eltype(A))

    other_sizes = ntuple(i -> i == dim_index ? 1 : size(A, i), ndims(A))

    for I in CartesianIndices(other_sizes)
        idx = ntuple(i -> i == dim_index ? Colon() : I[i], ndims(A))
        y = Array(A[idx...])
        out[idx...] = f(y)
    end

    return rebuild(A, out)
end

# -------------------------------
# Boxcar (moving average)
# -------------------------------
"""
    boxcar_filter_dim(A::AbstractDimArray, dim; window::Int)

Apply a boxcar (moving average) filter along the specified dimension.

# Arguments
- A: Input DimArray
- dim: Dimension along which to apply the filter
- window: Window size (must be odd)

# Returns
- Smoothed DimArray
"""
function boxcar_filter_dim(A::AbstractDimArray, dim; window::Int)
    @assert isodd(window) "window size must be odd"
    half = div(window, 2)

    function boxcar_1d(y)
        n = length(y)
        out = similar(y)
        for i = 1:n
            lo = max(1, i - half)
            hi = min(n, i + half)
            out[i] = mean(@view y[lo:hi])
        end
        return out
    end

    return _apply_along_dim(A, dim, boxcar_1d)
end

# -------------------------------
# Savitzky-Golay filter
# -------------------------------
"""
    sg_filter_dim(A::AbstractDimArray, dim;
                  window::Int, polyorder::Int)

Apply a Savitzky-Golay filter along the specified dimension.

# Arguments
- A: Input DimArray
- dim: Dimension along which to apply the filter
- window: Window size (must be odd)
- polyorder: Polynomial order (polyorder < window)

# Returns
- Smoothed DimArray
"""
function sg_filter_dim(A::AbstractDimArray, dim; window::Int, polyorder::Int)

    @assert isodd(window) "window size must be odd"
    @assert polyorder < window "polyorder must be < window"

    half = div(window, 2)

    # Precompute convolution coefficients
    x = collect((-half):half)
    X = hcat([x .^ k for k = 0:polyorder]...)
    # pseudo-inverse
    coeff = (X'X) \ X'
    # smoothing kernel (0th derivative at center)
    kernel = coeff[1, :]

    function sg_1d(y)
        n = length(y)
        out = similar(y)
        for i = 1:n
            lo = max(1, i - half)
            hi = min(n, i + half)

            # edge handling: shrink window
            xx = collect(lo:hi) .- i
            Xloc = hcat([xx .^ k for k = 0:polyorder]...)
            coeffloc = (Xloc'Xloc) \ Xloc'
            k = coeffloc[1, :]

            out[i] = sum(k .* y[lo:hi])
        end
        return out
    end

    return _apply_along_dim(A, dim, sg_1d)
end

# -------------------------------
# Binomial filter
# -------------------------------
"""
    binomial_filter_dim(A::AbstractDimArray, dim; order::Int)

Apply a binomial (Pascal triangle) smoothing filter along the specified dimension.

Kernel corresponds to coefficients of (1/2^n) * [C(n,0), ..., C(n,n)].

# Arguments
- A: Input DimArray
- dim: Dimension along which to apply the filter
- order: Binomial order (kernel size = order + 1)

# Returns
- Smoothed DimArray
"""
function binomial_filter_dim(A::AbstractDimArray, dim; order::Int)

    @assert order ≥ 1 "order must be ≥ 1"

    # Binomial coefficients
    coeffs = [binomial(order, k) for k = 0:order]
    kernel = coeffs ./ sum(coeffs)

    half = div(length(kernel), 2)

    function binomial_1d(y)
        n = length(y)
        out = similar(y)
        for i = 1:n
            acc = zero(eltype(y))
            wsum = zero(eltype(y))

            for j = 1:length(kernel)
                idx = i + (j - half - 1)
                if 1 ≤ idx ≤ n
                    w = kernel[j]
                    acc += w * y[idx]
                    wsum += w
                end
            end

            out[i] = acc / wsum
        end
        return out
    end

    return _apply_along_dim(A, dim, binomial_1d)
end

# -------------------------------
# Gaussian filter
# -------------------------------
"""
    gaussian_filter_dim(A::AbstractDimArray, dim; sigma::Real, radius::Int=ceil(Int, 3sigma))

Apply a Gaussian smoothing filter along the specified dimension.

# Arguments
- A: Input DimArray
- dim: Dimension along which to apply the filter
- sigma: Standard deviation of the Gaussian
- radius: Kernel radius (default ≈ 3σ)

# Returns
- Smoothed DimArray
"""
function gaussian_filter_dim(
    A::AbstractDimArray,
    dim;
    sigma::Real,
    radius::Int = ceil(Int, 3sigma),
)
    @assert sigma > 0 "sigma must be positive"

    x = collect((-radius):radius)
    kernel = exp.(-(x .^ 2) ./ (2sigma^2))
    kernel ./= sum(kernel)

    function gaussian_1d(y)
        n = length(y)
        out = similar(y)
        for i = 1:n
            acc = zero(eltype(y))
            wsum = zero(eltype(y))
            for j = 1:length(kernel)
                idx = i + (j - radius - 1)
                if 1 ≤ idx ≤ n
                    w = kernel[j]
                    acc += w * y[idx]
                    wsum += w
                end
            end
            out[i] = acc / wsum
        end
        return out
    end

    return _apply_along_dim(A, dim, gaussian_1d)
end

