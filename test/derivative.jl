using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using ARPES: _vector_diff, _gradient_modulus
using DimensionalData
using DimensionalData: Dim, hasdim, lookup
using Random
using Statistics

if !@isdefined(test_ARPESData)
    include("fixture.jl")
end

if !@isdefined(_random_monotonic_vector)
    function _random_monotonic_vector(n::Int, start = 0.0, stop = 10.0)
        @assert n > 1 "n must be greater than 1)"
        v = rand(n - 2)

        while length(unique(v)) < n - 2
            push!(v, rand())
            v = unique(v)
        end

        v = sort(vcat(0.0, 1.0, v))
        return v .* (stop - start) .+ start
    end
end

if !@isdefined(test_data)
    Random.seed!(123)

    n = 200
    t = range(0, 10, length = n)

    noise_level = 0.1
    signal_sine = sin.(t)
    noisy_sine = signal_sine .+ noise_level .* randn(n)

    gaussian(x, μ, σ) = exp.(-((x .- μ) .^ 2) ./ (2σ^2))
    signal_gauss = gaussian(t, 5.0, 0.5)
    noisy_gauss = signal_gauss .+ noise_level .* randn(n)

    test_data = (
        sine = (clean = signal_sine, noisy = noisy_sine),
        gauss = (clean = signal_gauss, noisy = noisy_gauss),
        t = t,
    )
end

@testset "Test for derivative 1st order" begin
    t_rnd = _random_monotonic_vector(200)
    Z_rnd = DimArray(sin.(t_rnd), (Dim{:t}(t_rnd)))
    Z = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t)))
    #
    Z_rnd_diff_1 = derivative(Z_rnd, :t)
    Z_diff_1 = derivative(Z, :t)
    cos_rnd = cos.(t_rnd)
    cos_clean = cos.(test_data.t)
    @test maximum(parent(Z_rnd_diff_1) - cos_rnd) < 0.005
    @test maximum(parent(Z_diff_1) - cos_clean) < 0.005
end

@testset "Test for _vector_diff and _gradient_modulus" begin
    test_array = [
        0 1 2 3 4 5;
        6 7 8 9 10 11;
        12 13 -2.5 15 16 17;
        18 19 20 21 2.5 23;
        24 25 26 27 28 29
    ]
    @test _vector_diff(test_array, (1, 0)) == [
        6.0 6.0 6.0 6.0 6.0 6.0;
        6.0 6.0 -10.5 6.0 6.0 6.0;
        6.0 6.0 22.5 6.0 -13.5 6.0;
        6.0 6.0 6.0 6.0 25.5 6.0
    ]
    @test isapprox(
        _vector_diff(test_array, (0, 1)),
        [
            1.0 1.0 1.0 1.0 1.0;
            1.0 1.0 1.0 1.0 1.0;
            1.0 -15.5 17.5 1.0 1.0;
            1.0 1.0 1.0 -18.5 20.5;
            1.0 1.0 1.0 1.0 1.0
        ],
    )
    @test isapprox(
        _vector_diff(test_array, (1, 1)),
        [
            7.0 7.0 7.0 7.0 7.0;
            7.0 -9.5 7.0 7.0 7.0;
            7.0 7.0 23.5 -12.5 7.0;
            7.0 7.0 7.0 7.0 26.5
        ],
    )
    @test isapprox(
        _vector_diff(test_array, (-1, -1)),
        [
            -7.0 -7.0 -7.0 -7.0 -7.0;
            -7.0 9.5 -7.0 -7.0 -7.0;
            -7.0 -7.0 -23.5 12.5 -7.0;
            -7.0 -7.0 -7.0 -7.0 -26.5;
        ],
    )
    @test isapprox(
        _vector_diff(test_array, (-1, 0)),
        [
            -6.0 -6.0 -6.0 -6.0 -6.0 -6.0;
            -6.0 -6.0 10.5 -6.0 -6.0 -6.0;
            -6.0 -6.0 -22.5 -6.0 13.5 -6.0;
            -6.0 -6.0 -6.0 -6.0 -25.5 -6.0;
        ],
    )

    @test isapprox(
        _gradient_modulus(test_array),
        [
            9.27362 10.583 10.583 10.583 10.583 7.87401;
            12.1244 16.225 17.2119 18.1452 14.8997 12.1244;
            12.1244 21.4767 48.9898 25.1893 19.1898 18.2277;
            12.1244 25.6759 26.3106 32.6573 57.1314 23.796;
            7.87401 10.583 10.583 26.2155 26.949 27.1892
        ];
        rtol = 1e-5,
    )
end
