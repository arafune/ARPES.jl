using Test
using DimensionalData
using Interfaces
using DimensionalData.Interfaces
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, CPS, TypeI, FinalStateEnergy
using ARPES.KConversion
using ARPES.KConversion: _deg2rad

@testset "Test for k_conversion with standard SPD data (minimum)" begin
    data = test_spd_standard()
    kx_range = -1.0:0.1:1.0
    ky_range = -1.0:0.1:1.0
    eV_range = 5.2:0.5:8.0
    converted_data = k_conversion(data)
    #    converted_data = k_conversion(
    #        data;
    #        kx_range = kx_range,
    #        #ky_range = ky_range,
    #        eV_range = eV_range,
    #        β0 = 10.0,
    #        ξ0 = 5.0,
    #        δ0 = 2.0,
    #        χ0 = 3.0,
    #    )
    @test converted_data isa ARPESData
    @test hasdim(converted_data, :kx)
    @test hasdim(converted_data, :eV)
    @test metadata(dims(converted_data, :eV))[:energy_definition] == FinalStateEnergy
end


