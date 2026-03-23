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
