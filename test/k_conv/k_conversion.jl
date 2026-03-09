using Test
using DimensionalData
using Interfaces
using DimensionalData.Interfaces
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, CPS, TypeI, FinalStateEnergy
using ARPES.KConversion
using ARPES.KConversion: _deg2rad, _check_arpesdata

if !@isdefined(test_spd_standard)
    include("fixture.jl")
end

@testset "Test for helper _deg2rad function" begin
    @test _deg2rad(0.0) == 0.0
    @test _deg2rad(180.0) ≈ π
    @test _deg2rad(90.0) ≈ π/2
    @test _deg2rad(0.0:2.0:180.0) ≈ 0.0:_deg2rad(2.0):_deg2rad(180.0)
    @test _deg2rad([-90.0, 0.0, 90.0]) ≈ [-π/2, 0, π/2]
end

@testset "Test for helper _check_arpesdata" begin
    data = test_spd_standard()
    @test_logs (:warn, "Missing metadata: χ. Using default value of 0.0 for conversion.") begin
        @test _check_arpesdata(data) == true
    end
    delete!(metadata(data), :workfunction)
    @test_throws ArgumentError _check_arpesdata(data)

    data = test_spd_standard()
    delete!(metadata(data), :β)
    @assert !haskey(metadata(data), :β)
    @assert !hasdim(data, :psi)
    @test_throws ArgumentError _check_arpesdata(data)

    data = test_spd_standard()
    old = dims(data)[dimnum(data, :eV)]
    new_dim = Dim{:energy}(lookup(old))
    data = rebuild(data, dims = Base.setindex(dims(data), new_dim, dimnum(data, :eV)))
    @test_throws ArgumentError _check_arpesdata(data)
    # Further tests can be added here to check the correctness of the conversion.
end
