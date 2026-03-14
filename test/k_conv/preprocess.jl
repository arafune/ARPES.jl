using Test
using DimensionalData
using Interfaces
using DimensionalData.Interfaces
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, CPS, TypeI, FinalStateEnergy
using ARPES.KConversion
using ARPES.KConversion: _kx_range, _ky_range
using ARPES.KConversion: _check_arpesdata, reshape_for_nd


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

@testset "TEst for reshape_for_nd" begin
    x = [1, 2, 3]
    y = reshape_for_nd(x, 2)
    @test size(y[1]) == (3,)
    @test y[1] == x

    x = 0.1:0.5:3
    y = 10:10:40
    z = -8:1.3:1
    result = reshape_for_nd(x, y, z)
    @test (length(result[1]), length(result[2]), length(result[3])) ==
          (length(x), length(y), length(z))
    @test result[1] == reshape(x, :, 1, 1)
    #@test result[1, 1, 1] == (x[1], y[1], z[1])
    #@test result[end, end, end] == (x[end], y[end], z[end])


end

