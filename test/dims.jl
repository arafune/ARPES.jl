using Test
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using DimensionalData

@testset "test energy_definition in eV metadata" begin
    @test metadata(dims(spd_standard)[dimnum(spd_standard, :eV)])[:energy_definition] ==
          FinalStateEnergy
end
