using Test
using DimensionalData
using ARPES

file_path = joinpath(pkgdir(ARPES), "testdata", "arpes_ch_resolved.itx")
metadata = Dict(:beta=>0.0)
arpes_ch_resolved = load(file_path, loc = "SPD", metadata = metadata)

@testset "ARPESData delegate" begin
    @test size(arpes_ch_resolved) == (40, 341, 24)
    @test dims(arpes_ch_resolved) ==
          (phi(-10.6875:0.5480769230769231:10.6875), eV(4.5:0.005:6.2), detector_ch(1:1:24))
    @test name(arpes_ch_resolved) == :Region1_1
    @test axes(arpes_ch_resolved) == (1:40, 1:341, 1:24)
    @test metadata(arpes_ch_resolved)[:lens_voltage] == "40V"
    @test arpes_ch_resolved[1, 1, 1] ≈ 67.8969
    @test arpes_ch_resolved.unit == CPS()
    @test metadata(arpes_ch_resolved)[:beta] == 0.0
end
