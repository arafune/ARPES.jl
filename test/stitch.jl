using Test
using DimensionalData
using ARPES
using ARPES: phi, eV, delay


function build_arrays()
    # Array A: X from 1 to 8
    x1 = range(1.0, stop = 8.0, length = 20)
    y1 = collect(range(5, stop = 9, length = 10))
    data1 = collect(1.0:200) |> d -> reshape(d, 20, 10)
    A = DimArray(data1, (X = x1, Y = y1))

    # Array B: X from 3 to 15 (Overlaps A)
    x2 = range(3.0, stop = 15.0, length = 20)
    y2 = collect(range(5, stop = 9, length = 10))
    data2 = collect(1.0:200) |> d -> reshape(d, 10, 20)
    B = DimArray(data2', (X = x2, Y = y2))

    # Array C: X from 13 to 15 (No overlap with A)
    x3 = range(13.0, stop = 15.0, length = 10)
    y3 = collect(range(5, stop = 9, length = 10))
    data3 = collect(1.0:100) |> d -> reshape(d, 10, 10)
    C = DimArray(data3', (X = x3, Y = y3))

    return A, B, C
end


function build_arrays3D()
    x = range(7.0, stop = 15.0, length = 20)
    y = collect(range(5, stop = 9, length = 10))
    z = collect(range(1, stop = 5, length = 5))
    data = collect(1.0:1000) * 2.1 |> d->reshape(d, 20, 10, 5)
    A = ARPESData(data, (phi = x, eV = y, delay = z))

    x = range(1.0, stop = 8.0, length = 20)
    y = collect(range(5, stop = 9, length = 10))
    z = collect(range(1, stop = 5, length = 5))
    data = collect(1.0:1000) |> d->reshape(d, 20, 10, 5)
    B = ARPESData(data, (phi = x, eV = y, delay = z))

    return A, B
end

@testset "ARPES Stitching Tests" begin
    A, B, C = build_arrays()

    @testset "Case: Overlapping Stitched (A & B)" begin
        res = stitch_along(A, B, :X, 0.5, 1.0)
        @test res isa DimArray
        @test size(res, :X) < (size(A, :X) + size(B, :X))
        @test minimum(lookup(res, :X)) == 1.0
        @test maximum(lookup(res, :X)) == 15.0
        @test issorted(lookup(res, :X))
    end

    @testset "Case: Disjoint Stitched (A & C)" begin
        res = stitch_along(A, C, :X, 0.5, 1.0)
        @test size(res, :X) == size(A, :X) + size(C, :X)
        @test issorted(lookup(res, :X))
    end

    @testset "Case: Intensity Enhancement" begin
        gain = 2.5
        res = stitch_along(A, C, :X, nothing, gain)
        @test res[X(At(1.0)), Y(At(5.0))] == A[X(At(1.0)), Y(At(5.0))] * gain
        @test res[X(At(15.0)), Y(At(5.0))] == C[X(At(15.0)), Y(At(5.0))]
    end

    @testset "Case: Error Handling" begin
        @test_throws AssertionError stitch_along(A, B, :X, 1.2, 1.0)

        B_wrong = B[Y(1:5)]
        @test_throws AssertionError stitch_along(A, B_wrong, :X, 0.5, 1.0)

        @test_throws ErrorException stitch_along(A, A, :X, 0.5, 1.0)
    end
end

@testset "ARPES Stitching Tests - 3D" begin
    A, B = build_arrays3D()

    @testset "Case: Overlapping Stitched (A & B)" begin
        res = stitch_along(A, B, :phi; seam_ratio = 0.5, gain_a = 1.0)
        @test res isa ARPESData
        @test size(res, :phi) < (size(A, :phi) + size(B, :phi))
        @test minimum(lookup(res, :phi)) == 1.0
        @test maximum(lookup(res, :phi)) == 15.0
        @test issorted(lookup(res, :phi))
    end
end
