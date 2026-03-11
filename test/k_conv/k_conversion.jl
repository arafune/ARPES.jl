using Test
using DimensionalData
using Interfaces
using DimensionalData.Interfaces
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, CPS, TypeI, FinalStateEnergy
using ARPES.KConversion
using ARPES.KConversion: _deg2rad

if !@isdefined(test_spd_standard)
    include("fixture.jl")
end

@testset "Test for helper _deg2rad function" begin
    @test _deg2rad(0.0) == 0.0
    @test _deg2rad(180.0) ≈ π
    @test _deg2rad(90.0) ≈ π/2
    r = _deg2rad(0.0:2.0:180.0, 1.0)
    @test first(r) == _deg2rad(0.0, 1.0)
    @test step(r) ≈ _deg2rad(2.0, 0.0)
    @test last(r) ≈ _deg2rad(180.0, 1.0)
    @test _deg2rad([-90.0, 0.0, 90.0]) ≈ [-π/2, 0, π/2]
end

