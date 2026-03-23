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
        P = Matrix{Float64}(LinearAlgebra.I, nx, nx) * 1e3

        x_est = zeros(nx, n)

        for t = 1:n
            dt = Δt[t]

            if model == :rw
                F = reshape([1.0], 1, 1)
                Q = reshape([q * max(dt, 1e-12)], 1, 1)
                H = reshape([1.0], 1, 1)
            else # :cv
                F = [
                    1.0 dt;
                    0.0 1.0
                ]
                Q = q * [
                    dt^3/3 dt^2/2;
                    dt^2/2 dt
                ]
                H = reshape([1.0 0.0], 1, 2)
            end

            # Predict
            x̂ = F * x̂
            P = F * P * F' + Q

            # Update (skip NaN)
            if !isnan(y[t])
                S = (H*P*H')[1] + R
                K = P * H' / S
                residual = y[t] - (H*x̂)[1]
                x̂ = x̂ + K * residual
                P = (LinearAlgebra.I - K * H) * P
            end

            x_est[:, t] = x̂
        end

        # Write back (use first state component)
        out[idx...] = x_est[1, :]
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
function boxcar_filter_dim(A::AbstractDimArray, dim; window::Real, pixel = true)
    if pixel
        n_pixels = Int(window)
    else
        s = step(dims(A, dim))
        n_pixels = round(Int, window / s)
        n_pixels = iseven(n_pixels) ? n_pixels + 1 : n_pixels
    end
    @assert isodd(n_pixels) "Compute pixel window size ($n_pixels) must be odd"
    @assert n_pixels >= 1 "Window size is too small for the given dimension step"

    half = div(n_pixels, 2)

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
function sg_filter_dim(
    A::AbstractDimArray,
    dim;
    window::Real,
    polyorder::Int = 2,
    pixel = true,
)
    # 1. Determine n_pixels (must be odd)
    if pixel
        n_pixels = Int(window)
    else
        s = abs(step(dims(A, dim)))
        n_pixels = round(Int, window / s)
        n_pixels = iseven(n_pixels) ? n_pixels + 1 : n_pixels
    end

    @assert isodd(n_pixels) "Window size ($n_pixels) must be odd"
    @assert polyorder < n_pixels "polyorder ($polyorder) must be less than window ($n_pixels)"

    half = div(n_pixels, 2)

    # 2. Precompute the Vandermonde matrix and its pseudo-inverse for a fixed window
    x_grid = collect((-half):half)
    X = hcat([x_grid .^ k for k = 0:polyorder]...)
    # The pseudo-inverse matrix (polyorder+1 x n_pixels)
    # Each row `m` gives the coefficients for the `(m-1)`-th derivative
    B = (X'X) \ X'

    function sg_1d(y)
        n = length(y)
        out = similar(y)
        for i = 1:n
            # Find the window that includes index i. 
            # Clamp the center so the window [center-half : center+half] stays within 1:n.
            center = clamp(i, half + 1, n - half)
            lo, hi = center - half, center + half

            # Distance from the window center to the current point i
            # (ranges from -half at the very left edge to +half at the right edge)
            relative_pos = i - center

            # --- The Key Fix ---
            # To get the smoothed value at `relative_pos`, we need to evaluate 
            # the fitted polynomial: P(x) = c0 + c1*x + c2*x^2 + ...
            # This is equivalent to taking a weighted sum of the window, 
            # where weights are: w = B[1,:] + B[2,:]*x + B[3,:]*x^2 + ...

            weights = zeros(eltype(B), n_pixels)
            for k = 0:polyorder
                # B[k+1, :] gives the weights for the k-th derivative term
                # We scale it by (relative_pos^k / k!) if we were doing derivatives, 
                # but for polynomial evaluation at `relative_pos`:
                weights .+= B[k+1, :] .* (relative_pos^k)
            end

            out[i] = sum(weights .* y[lo:hi])
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
    gaussian_filter_dim(A::AbstractDimArray, dim; sigma::Real, radius::Int=ceil(Int, 3sigma), pixel=true)

Apply a Gaussian smoothing filter along the specified dimension.

If `pixel=false`, `sigma` is interpreted in physical units of the dimension.

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
    radius::Int = 0,
    pixel = true,
)
    @assert sigma > 0 "sigma must be positive"

    # 1. Convert physical sigma to pixel-based sigma if necessary
    sigma_pix = if pixel
        Float64(sigma)
    else
        # Get the absolute step size of the dimension (e.g., eV or degrees)
        s = abs(step(dims(A, dim)))
        sigma / s
    end

    # 2. Determine the kernel radius
    # Default to 3 * sigma to cover 99.7% of the Gaussian distribution
    r = radius > 0 ? radius : ceil(Int, 3 * sigma_pix)

    # 3. Create the Gaussian kernel
    x = collect((-r):r)
    kernel = exp.(-(x .^ 2) ./ (2 * sigma_pix^2))
    # Pre-normalize, though we will re-normalize at boundaries
    kernel ./= sum(kernel)

    # Internal function to process a 1D slice
    function gaussian_1d(y)
        n = length(y)
        out = similar(y)
        for i = 1:n
            acc = zero(eltype(y))
            wsum = zero(eltype(y))
            for j = 1:length(kernel)
                # Data index corresponding to the kernel index j
                idx = i + x[j]

                if 1 ≤ idx ≤ n
                    w = kernel[j]
                    acc += w * y[idx]
                    wsum += w
                end
            end
            # Re-normalize to handle boundary truncation correctly
            out[i] = acc / wsum
        end
        return out
    end

    return _apply_along_dim(A, dim, gaussian_1d)
end

