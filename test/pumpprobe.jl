using Test
using DimensionalData
using DimensionalData: @dim, Dim
using ARPES
using ARPES: ARPESData, kx, ky, kz, phi, psi, eV, delay

@testset "ARPESPlots tarpes Heatmap Tests" begin
    data = rand(Float64, 40, 60, 30)
    data[3, 3, 3] = NaN
    data[1, 1, 1] = NaN
    A = ARPESData(
        data,
        (phi(range(-10, 10.0, 40)), eV(range(0, 5.0, 60)), delay(range(-3.0, 5.1, 30))),
    )

    # Basic scalar selection: snapshot shape and temporal evolution size
    arpes_snapshot, temporal = tarpes_evolution(A, 0.0, 0.0)
    @test size(arpes_snapshot) == (size(A, 1), size(A, 2))

    delays = dims(A, :delay)
    expected_rows = sum(delays .<= 0.0)
    @test size(temporal, 1) == expected_rows
    @test size(temporal, 2) == size(arpes_snapshot, 2)

    # delay-index variant should match the delay_time variant at the nearest index
    nearest_idx = argmin(abs.(dims(A, :delay) .- 0.0))
    arpes_snapshot_idx, temporal_idx = tarpes_evolution(A, nearest_idx, 0.0)
    @test isequal(arpes_snapshot, arpes_snapshot_idx)
    @test isequal(temporal, temporal_idx)

    # Tuple averaging along the non-dispersion axis
    center = 0.0
    width = 2.0
    _, temporal_avg = tarpes_evolution(A, 0.0, (center, width))

    # Build expected temporal evolution manually
    expected =
        A[Dim{:phi}(Between(center - width/2, center + width/2))] |>
        x ->
            dropdims(x, dims = :phi) |>
            x ->
                mean(x, dims = :phi) |>
                x -> permutedims(x, (:delay, :eV)) |> x[Dim{:delay}(Between(-Inf, 0.0))]

    @test size(temporal_avg) == size(expected)
    @test isapprox(temporal_avg, expected; atol = 1e-12, rtol = 1e-8)

end
