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
    @test _check_arpesdata(data) == true
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

@testset "Test for reshape_for_nd" begin
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
end

@testset "Test _kx_range _ky_range" begin
    ek = 0.7990000000000004:0.002:1.9990000000000003
    α = 0.19966566642815128:-0.0020943951023931952:-0.21711895894809455
    β_ = 0.0
    χ_ = 0.0
    ξ_ = -0.0
    δ_ = 0.0

    kx_range = _kx_range(TypeI, α, β_, ek, χ_, ξ_, δ_)
    @test kx_range[1] ≈ -0.1436676048385842
    @test step(kx_range) ≈ 0.0009368123759122965
    @test kx_range[end] ≈ 0.15517554307743836
    @test _ky_range(TypeI, α, β_, ek, χ_, ξ_, δ_) == 0.0

end
