using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using ARPES: _apply_along_dim
using DimensionalData
using DimensionalData: Dim, hasdim, lookup
using Random
using Statistics

# 1. Basic setting
#
#using Random

Random.seed!(123)

n = 200
t = range(0, 10, length = n)
noise_level = 0.1

# -------------------------------
# Signals
# -------------------------------

# Sine wave
signal_sine = sin.(t)
noisy_sine = signal_sine .+ noise_level .* randn(n)

# Gaussian pulse
gaussian(x, μ, σ) = exp.(-((x .- μ) .^ 2) ./ (2σ^2))

signal_gauss = gaussian(t, 5.0, 0.5)
noisy_gauss = signal_gauss .+ noise_level .* randn(n)

# -------------------------------
# Test data
# -------------------------------

test_data = (
    sine = (clean = signal_sine, noisy = noisy_sine),
    gauss = (clean = signal_gauss, noisy = noisy_gauss),
    t = t,
)

@testset "boxcar_filter_dim" begin
    A = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    B = boxcar_filter_dim(A, :t, window = 5)

    @test size(B) == size(A)
    @test dims(B) == dims(A)

    @test var(B) < var(A)
end

@testset "gaussian_filter_dim" begin
    A = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    B = gaussian_filter_dim(A, :t, sigma = 2.0)

    @test size(B) == size(A)

    # smoothness check
    @test var(diff(parent(B))) < var(diff(parent(A)))
end

@testset "binomial_filter_dim" begin
    A = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    B = binomial_filter_dim(A, :t, order = 4)

    @test size(B) == size(A)
    @test var(B) < var(A)
end

@testset "sg_filter_dim polynomial preservation" begin
    t = range(0, 10, length = 100)
    y = 3 .+ 2t .+ 0.5t .^ 2

    A = DimArray(y, (Dim{:t}(t),))
    B = sg_filter_dim(A, :t, window = 7, polyorder = 2)

    @test parent(B) ≈ parent(A) atol=1e-10
end

@testset "kalman_smooth_dim_llpf" begin
    using Random
    Random.seed!(123)

    noise = 0.1 .* randn(length(t))
    signal = sin.(t)
    noisy = signal .+ noise

    A = DimArray(noisy, (Dim{:t}(t),))

    B = kalman_smooth_dim_llpf(A, :t, R = 0.01, q = 0.1)

    # shape
    @test size(B) == size(A)

    # smoothing effect
    @test var(parent(B) .- signal) < var(parent(A) .- signal)
end

@testset "kalman non-uniform dt" begin
    t = sort(rand(100) .* 10)
    signal = sin.(t)

    A = DimArray(signal, (Dim{:t}(t),))
    B = kalman_smooth_dim_llpf(A, :t)

    @test size(B) == size(A)
end

@testset "nan handling" begin
    t = sort(rand(100) .* 10)
    signal = sin.(t)
    y = copy(signal)
    y[50] = NaN

    A = DimArray(y, (Dim{:t}(t),))
    B = kalman_smooth_dim_llpf(A, :t)

    @test !all(isnan, parent(B))
end

@testset "_apply_along_dim" begin

    # -------------------------------
    # 1. 1D: identity (no-op) test
    # Ensures data and metadata are preserved
    # -------------------------------
    t = range(0, 10, length = 50)
    A = DimArray(collect(t), (Dim{:t}(t),))

    B = _apply_along_dim(A, :t, x -> x)

    @test parent(B) == parent(A)
    @test dims(B) == dims(A)
    @test size(B) == size(A)

    # -------------------------------
    # 2. 1D: simple transformation
    # Verify function is correctly applied
    # -------------------------------
    B = _apply_along_dim(A, :t, x -> x .+ 1)

    @test parent(B) ≈ parent(A) .+ 1

    # -------------------------------
    # 3. 2D: apply along first dimension
    # Each column should be processed independently
    # -------------------------------
    x = 1:4
    y = 1:3
    data = reshape(collect(1.0:12.0), 4, 3)

    A2 = DimArray(data, (Dim{:x}(x), Dim{:y}(y)))

    # Subtract mean along :x (column-wise)
    B2 = _apply_along_dim(A2, :x, v -> v .- mean(v))

    for j in axes(data, 2)
        @test mean(parent(B2)[:, j]) ≈ 0
    end

    @test size(B2) == size(A2)
    @test dims(B2) == dims(A2)

    # -------------------------------
    # 4. Apply along second dimension
    # Each row should be processed independently
    # -------------------------------
    B3 = _apply_along_dim(A2, :y, v -> v .- mean(v))

    for i in axes(data, 1)
        @test mean(parent(B3)[i, :]) ≈ 0
    end

    # -------------------------------
    # 5. Position-dependent function
    # Ensures correct indexing inside each slice
    # -------------------------------
    B4 = _apply_along_dim(A2, :x, v -> v .* (1:length(v)))

    for j in axes(data, 2)
        expected = data[:, j] .* (1:length(x))
        @test parent(B4)[:, j] == expected
    end

    # -------------------------------
    # 6. Type stability
    # Output element type should be preserved
    # -------------------------------
    Afloat = DimArray(rand(Float32, 5, 5), (Dim{:x}(1:5), Dim{:y}(1:5)))
    Bfloat = _apply_along_dim(Afloat, :x, x -> x .* 2)

    @test eltype(Bfloat) == Float32

    # -------------------------------
    # 7. NaN propagation
    # Function should not silently remove NaNs
    # -------------------------------
    A_nan = DimArray([1.0, NaN, 3.0], (Dim{:t}(1:3),))
    B_nan = _apply_along_dim(A_nan, :t, x -> x)

    @test isnan(parent(B_nan)[2])

    # -------------------------------
    # 8. 3D test
    # Ensure correct behavior in higher dimensions
    # -------------------------------
    A3 = DimArray(rand(4, 3, 2), (Dim{:x}(1:4), Dim{:y}(1:3), Dim{:z}(1:2)))

    B3d = _apply_along_dim(A3, :x, v -> v .- mean(v))

    for j = 1:3, k = 1:2
        @test mean(parent(B3d)[:, j, k]) ≈ 0 atol=1e-12
    end

    @test size(B3d) == size(A3)
    @test dims(B3d) == dims(A3)

    # -------------------------------
    # 9. Consistency with mapslices
    # Compare with Base.mapslices behavior
    # -------------------------------
    f(v) = v .- mean(v)

    B_ref = mapslices(f, parent(A2); dims = 1)  # dims=1 corresponds to :x
    B_test = _apply_along_dim(A2, :x, f)

    @test parent(B_test) ≈ B_ref

end
