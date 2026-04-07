using Test
using DimensionalData
using Interfaces
using DimensionalData.Interfaces
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, CPS, TypeI, FinalStateEnergy
using ARPES.KConversion
using ARPES.KConversion: _deg2rad

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


#@teset "Test for k_conversion with standard SPD data (with parameters)" begin
#    converted_data_2 = k_conversion(data; β0 = 10.0, ξ0 = 5.0, δ0 = 2.0, χ0 = 3.0)
#    @test converted_data_2 isa ARPESData
#    @test hasdim(converted_data_2, :kx)
#    @test hasdim(converted_data_2, :eV)
#    @test metadata(dims(converted_data_2, :eV))[:energy_definition] == FinalStateEnergy
#end
