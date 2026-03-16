using Test
using ARPES.KConversion: _interpolate
@testset "_interpolate tests" begin
    # 2D tests
    c0 = 1.0:1.0:5.0
    c1 = 10.0:10.0:50.0
    values = [i + j for i in c0, j in c1]  # 5x5 matrix

    p0 = [2.5, 4.0]
    p1 = [15.0, 35.0]
    res = _interpolate(p0, p1, c0, c1, values)
    @test length(res) == 2
    @test isapprox(res[1], 2.5 + 15.0; atol = 1e-8)
    @test isapprox(res[2], 4.0 + 35.0; atol = 1e-8)

    c0_array = collect(c0)
    c1_array = collect(c1)


    res1_1 = _interpolate(p0, p1, c0_array, c1_array, values)
    @test length(res1_1) == 2
    @test isapprox(res, res1_1, atol = 1e-8)

    res2 = _interpolate([0.0], [0.0], c0, c1, values)
    @test isnan(res2[1])

    c0r = 5.0:-1.0:1.0
    valuesr = reverse(values, dims = 1)
    res3 = _interpolate([2.0], [20.0], c0r, c1, valuesr)
    @test isapprox(res3[1], 2.0 + 20.0; atol = 1e-8)

    c1r = 50.0:-10.0:10.0
    valuesr2 = reverse(values, dims = 2)
    res4 = _interpolate([3.0], [30.0], c0, c1r, valuesr2)
    @test isapprox(res4[1], 3.0 + 30.0; atol = 1e-8)

    valuesr3 = reverse(reverse(values, dims = 1), dims = 2)
    res5 = _interpolate([4.0], [40.0], c0r, c1r, valuesr3)
    @test isapprox(res5[1], 4.0 + 40.0; atol = 1e-8)

    res6 = _interpolate([0.0], [60.0], c0r, c1r, valuesr3)
    @test isnan(res6[1])

    res7 = _interpolate([1.0], [10.0], c0, c1, values)
    @test isapprox(res7[1], 1.0 + 10.0; atol = 1e-8)

    c0n = 1.0:1.0:3.0
    c1n = 10.0:10.0:40.0
    valuesn = [i + j for i in c0n, j in c1n]  # 3x4 matrix
    res8 = _interpolate([2.0], [20.0], c0n, c1n, valuesn)
    @test isapprox(res8[1], 2.0 + 20.0; atol = 1e-8)

    dummy1 = 0.3
    dummy2 = 0.4
    result_pseudo_3d_1 = _interpolate(p0, p1, dummy1, c0, c1, dummy2, values)
    result_pseudo_3d_2 = _interpolate(p0, dummy1, p1, c0, dummy2, c1, values)
    result_pseudo_3d_3 = _interpolate(dummy1, p0, p1, dummy2, c0, c1, values)
    @test isapprox(result_pseudo_3d_1, res, atol = 1e-8)
    @test isapprox(result_pseudo_3d_2, res, atol = 1e-8)
    @test isapprox(result_pseudo_3d_3, res, atol = 1e-8)

    result_pseudo_3d_4 = _interpolate(p0, p1, dummy1, c0_array, c1_array, dummy2, values)
    result_pseudo_3d_5 = _interpolate(p0, dummy1, p1, c0, dummy2, c1_array, values)
    result_pseudo_3d_6 = _interpolate(dummy1, p0, p1, dummy2, c0_array, c1_array, values)
    @test isapprox(result_pseudo_3d_4, res, atol = 1e-8)
    @test isapprox(result_pseudo_3d_5, res, atol = 1e-8)
    @test isapprox(result_pseudo_3d_6, res, atol = 1e-8)

    # 3D tests
    c0_3d = 1.0:1.0:4.0
    c1_3d = 10.0:10.0:40.0
    c2_3d = 100.0:100.0:400.0
    values3d = [i + j + k for i in c0_3d, j in c1_3d, k in c2_3d]  # 4x4x4 array

    p0_3d = [2.5, 3.0]
    p1_3d = [15.0, 25.0]
    p2_3d = [150.0, 350.0]
    res3d = _interpolate(p0_3d, p1_3d, p2_3d, c0_3d, c1_3d, c2_3d, values3d)
    @test length(res3d) == 2
    @test isapprox(res3d[1], 2.5 + 15.0 + 150.0; atol = 1e-8)
    @test isapprox(res3d[2], 3.0 + 25.0 + 350.0; atol = 1e-8)

    c0_3d_array = collect(c0_3d)
    c1_3d_array = collect(c1_3d)
    c2_3d_array = collect(c2_3d)
    res3d_1 =
        _interpolate(p0_3d, p1_3d, p2_3d, c0_3d_array, c1_3d_array, c2_3d_array, values3d)
    @test length(res3d_1) == 2
    @test isapprox(res3d, res3d_1, atol = 1e-8)


    # Out-of-bounds for 3D
    res3d_oob = _interpolate([0.0], [0.0], [0.0], c0_3d, c1_3d, c2_3d, values3d)
    @test isnan(res3d_oob[1])

    # Reversed axes for 3D
    c0r_3d = 4.0:-1.0:1.0
    c1r_3d = 40.0:-10.0:10.0
    c2r_3d = 400.0:-100.0:100.0
    values3d_r = reverse(reverse(reverse(values3d, dims = 1), dims = 2), dims = 3)
    res3d_rev = _interpolate([2.0], [20.0], [200.0], c0r_3d, c1r_3d, c2r_3d, values3d_r)
    @test isapprox(res3d_rev[1], 2.0 + 20.0 + 200.0; atol = 1e-8)

end
