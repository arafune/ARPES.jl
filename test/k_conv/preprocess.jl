using Test
using DimensionalData
using Interfaces
using DimensionalData.Interfaces
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, CPS, TypeI, FinalStateEnergy
using ARPES.KConversion
using ARPES.KConversion: _kx_range, _ky_range
using ARPES.KConversion: _check_arpesdata


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


