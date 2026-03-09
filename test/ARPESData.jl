using Test
using DimensionalData
using Interfaces
using DimensionalData.Interfaces
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, CPS, TypeI, FinalStateEnergy

if !@isdefined(test_spd_standard)
    include("fixture.jl")
end


@testset "ARPESData delegate" begin
    arpes_ch_resolved = test_arpes_ch_resolved()
    @test size(arpes_ch_resolved) == (40, 341, 24)
    @test dims(arpes_ch_resolved) ==
          (phi(-10.6875:0.5480769230769231:10.6875), eV(4.5:0.005:6.2), detector_ch(1:1:24))
    @test name(arpes_ch_resolved) == :Region1_1
    @test axes(arpes_ch_resolved) == (1:40, 1:341, 1:24)
    @test metadata(arpes_ch_resolved)[:lens_voltage] == "40V"
    @test arpes_ch_resolved[1, 1, 1] ≈ 67.8969
    @test metadata(arpes_ch_resolved)[:intensity_unit] == CPS
    @test metadata(arpes_ch_resolved)[:beta] == 0.0
    @test typeof(parent(arpes_ch_resolved)) == Array{Float64,3}

    @test eltype(arpes_ch_resolved) == Float64
    @test ndims(arpes_ch_resolved) == 3

    @test Base.IndexStyle(typeof(arpes_ch_resolved)) ==
          Base.IndexStyle(typeof(arpes_ch_resolved.data))

    b1 = Base.Broadcast.broadcastable(arpes_ch_resolved)
    b2 = Base.Broadcast.broadcastable(arpes_ch_resolved.data)
    @test typeof(b1) == typeof(b2)
    @test metadata(arpes_ch_resolved)[:hv] ≈ 4.835
end

@testset "test energy_definition in eV metadata" begin
    spd_standard = test_spd_standard()
    @test metadata(dims(spd_standard)[dimnum(spd_standard, :eV)])[:energy_definition] ==
          FinalStateEnergy
end

@testset "DimArrayInterface" begin
    spd_standard = test_spd_standard()
    Interfaces.test(DimensionalData.DimArrayInterface, spd_standard)
end
