using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using ARPES: _apply_along_dim
using DimensionalData
using DimensionalData: Dim, hasdim, lookup
using Random
using Statistics

# Basic setting

function _random_monotonic_vector(n::Int, start = 0.0, stop = 10.0)
    @assert n>1 "n must be greater than 1)"
    v = rand(n-2)

    while length(unique(v)) < n-2
        push!(v, rand())
        v = unique(v)
    end

    v = sort(vcat(0.0, 1.0, v))
    return v .* (stop - start) .+ start
end

Random.seed!(123)

n = 200
t = range(0, 10, length = n)
t_rnd = _random_monotonic_vector(n, 0, 10)

noise_level = 0.1

# Sine wave
signal_sine = sin.(t)
noisy_sine = signal_sine .+ noise_level .* randn(n)

# Gaussian pulse
gaussian(x, μ, σ) = exp.(-((x .- μ) .^ 2) ./ (2σ^2))

signal_gauss = gaussian(t, 5.0, 0.5)
noisy_gauss = signal_gauss .+ noise_level .* randn(n)

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

@testset "boxcar_filter_dim with noise" begin
    Z = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    A = DimArray(test_data.sine.noisy, (Dim{:t}(test_data.t),))
    B = boxcar_filter_dim(A, :t, window = 5)

    @test size(B) == size(A)
    @test dims(B) == dims(A)

    @test sum((Z-B) .^ 2) < sum((Z-A) .^ 2)
end

@testset "boxcar_filter_dim with noise" begin
    Z = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    A = DimArray(test_data.sine.noisy, (Dim{:t}(test_data.t),))
    B = boxcar_filter_dim(A, :t, window = 0.1, pixel = false)

    @test size(B) == size(A)
    @test dims(B) == dims(A)

    @test sum((Z-B) .^ 2) < sum((Z-A) .^ 2)
end



@testset "gaussian_filter_dim" begin
    A = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    B = gaussian_filter_dim(A, :t, sigma = 2.0)

    @test size(B) == size(A)

    # smoothness check
    @test var(diff(parent(B))) < var(diff(parent(A)))
end


@testset "gaussian_filter_dim with noise" begin
    Z = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    A = DimArray(test_data.sine.noisy, (Dim{:t}(test_data.t),))
    B = gaussian_filter_dim(A, :t, sigma = 5)

    @test size(B) == size(A)
    @test dims(B) == dims(A)

    @test sum((Z-B) .^ 2) < sum((Z-A) .^ 2)
end

@testset "gaussian_filter_dim with noise with physical value" begin
    Z = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    A = DimArray(test_data.sine.noisy, (Dim{:t}(test_data.t),))
    B = gaussian_filter_dim(A, :t, sigma = 0.05, pixel = false)

    @test size(B) == size(A)
    @test dims(B) == dims(A)

    @test sum((Z-B) .^ 2) < sum((Z-A) .^ 2)
end



@testset "binomial_filter_dim" begin
    A = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    B = binomial_filter_dim(A, :t, order = 4)

    @test size(B) == size(A)
    @test var(B) < var(A)
end

@testset "binomial_filter_dim with noise" begin
    Z = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    A = DimArray(test_data.sine.noisy, (Dim{:t}(test_data.t),))
    B = binomial_filter_dim(A, :t, order = 4)

    @test size(B) == size(A)
    @test dims(B) == dims(A)

    @test sum((Z-B) .^ 2) < sum((Z-A) .^ 2)
end

@testset "sg_filter_dim polynomial preservation" begin
    t = range(0, 10, length = 100)
    y = 3 .+ 2t .+ 0.5t .^ 2

    A = DimArray(y, (Dim{:t}(t),))
    B = sg_filter_dim(A, :t, window = 7, polyorder = 2)

    @test parent(B) ≈ parent(A) atol=1e-10
end

@testset "sg_filter_dim with noise" begin
    Z = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t),))
    A = DimArray(test_data.sine.noisy, (Dim{:t}(test_data.t),))
    B = sg_filter_dim(A, :t, window = 7, polyorder = 2)

    @test size(B) == size(A)
    @test dims(B) == dims(A)

    @test sum((Z-B) .^ 2) < sum((Z-A) .^ 2)
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

@testset "kalman_smooth_dim_llpf with noise" begin
    t_rnd = _random_monotonic_vector(n)
    Z_rnd = DimArray(sin.(t_rnd), (Dim{:t}(t_rnd)))
    A_rnd = DimArray(sin.(t_rnd) .+ noise_level .* randn(n), (Dim{:t}(t_rnd)))
    B_rnd = kalman_smooth_dim_llpf(A_rnd, :t, R = 0.01, q = 0.1)

    @test size(B_rnd) == size(A_rnd)
    @test sum((Z_rnd-B_rnd) .^ 2) < sum((Z_rnd-A_rnd) .^ 2)
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
