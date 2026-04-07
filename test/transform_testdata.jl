using Test
using ARPES
using ARPES: add_dim
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using ARPES: _shift_irregular
using DimensionalData
using DimensionalData: Dim, hasdim, lookup, Bins

if !@isdefined(test_spd_standard)
    include("fixture_testdata.jl")
end


@testset "Test for rebin" begin
    spd_standard = test_spd_standard()
    # Rebin along phi with new edges
    rebin_spd = rebin(spd_standard, :phi, 10)
    @test size(rebin_spd) == (10, 601)
    @test rebin_spd == rebin(spd_standard, :phi, Bins(10))
    @test ARPES._is_equal_spacing(lookup(rebin_spd, :phi))
end
