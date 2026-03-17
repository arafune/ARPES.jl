using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using DimensionalData: Dim, hasdim


@testset "Test for merge_consensus" begin
    d1 = Dict(:a => 1, :b => 2)
    d2 = Dict(:a => 1, :b => 3)
    @test merge_consensus(d1, d2) == Dict(:a => 1)
end
