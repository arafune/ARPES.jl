using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using DimensionalData: Dim, hasdim


@testset "Test for merge_consensus" begin
    d1 = Dict(:a => 1, :b => 2)
    d2 = Dict(:a => 1, :b => 3)
    @test merge_consensus(d1, d2) == Dict(:a => 1)
end

@testset "Test for merge_consensus with conflicted key" begin
    d1 = Dict(:a => 1, :b => 2)
    d2 = Dict(:a => 2, :b => 2)  # :a is different from d1
    d3 = Dict(:a => 3, :b => 2)  # :a appears again
    @test merge_consensus(d1, d2, d3) == Dict(:b => 2)
end

