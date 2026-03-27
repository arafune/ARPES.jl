using Test
using ARPES
using ARPES: add_dim, cat_arpes
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using ARPES: _shift_irregular
using DimensionalData
using DimensionalData: Dim, hasdim, lookup, Bins

if !@isdefined(test_ARPESData)
    include("fixture.jl")
end


@testset "Test for rebuild_with_slice" begin
    # --- setup ---
    x = range(1.0, stop = 3.0, length = 3)
    y = collect(range(5, stop = 8, length = 4))
    z = [10.0, 10.9, 12.0, 13.0, 14.0]

    data = reshape(collect(1.0:60.0), 3, 4, 5)
    A = DimArray(data, (X = x, Y = y, Z = z))

    # valid selector (fix Z dimension)
    dimsel = Z(At(z[3]))  # -> 3rd slice

    # -----------------------------
    # 1. Normal case (AbstractArray)
    # -----------------------------
    target = A[dimsel]
    X = fill(999.0, size(target))

    B = rebuild_with_slice(A, dimsel, X)

    @test B !== A                       # ensure non-mutating
    @test parent(B[dimsel]) == X        # values replaced correctly
    @test dims(B) == dims(A)            # metadata preserved

    # -----------------------------
    # 2. Normal case (AbstractDimArray)
    # -----------------------------
    X_dim = DimArray(X, dims(target))

    B2 = rebuild_with_slice(A, dimsel, X_dim)

    @test parent(B2[dimsel]) == X       # values replaced correctly
    @test dims(B2) == dims(A)           # metadata preserved

    # -----------------------------
    # 3. Size mismatch (both methods)
    # -----------------------------
    badX = fill(1.0, size(target, 1), size(target, 2) + 1)

    @test_throws DimensionMismatch rebuild_with_slice(A, dimsel, badX)
    @test_throws DimensionMismatch rebuild_with_slice(
        A,
        dimsel,
        DimArray(badX, dims(target)),
    )

    # -----------------------------
    # 4. Dimension metadata mismatch (DimArray only)
    # -----------------------------
    tdims = dims(target)

    # modify only coordinate values of the first dimension
    d1 = tdims[1]
    d2 = tdims[2]

    # extract raw values
    vals = collect(d1)

    # shift values slightly
    wrong_vals = vals .+ 1e-6

    # reconstruct dimension properly
    wrong_d1 = rebuild(d1, wrong_vals)

    wrong_dims = (wrong_d1, d2)

    X_bad_dim = DimArray(X, wrong_dims)

    @test_throws ArgumentError rebuild_with_slice(A, dimsel, X_bad_dim)

    # -----------------------------
    # 5. Invalid dimsel: fixes no dimension
    # -----------------------------
    bad_sel0 = (:)

    @test_throws ArgumentError rebuild_with_slice(A, bad_sel0, X)

    # -----------------------------
    # 6. Invalid dimsel: fixes multiple dimensions
    # -----------------------------
    bad_sel2 = (DimensionalData.X(At(x[1])), DimensionalData.Y(At(y[1])))

    @test_throws ArgumentError rebuild_with_slice(A, bad_sel2, X)

end


@testset "test for shift (1)" begin
    test_dim2D_1 = test_DimArray2D()
    test_dim3D_1 = test_DimArray3D()
    shift_2d_1 = shift(test_dim2D_1, :Y, -0.5)
    @test isequal(parent(shift_2d_1), [2.5 5.5 8.5 NaN; 3.5 6.5 9.5 NaN; 4.5 7.5 10.5 NaN])

    shift_2d_2 = _shift_irregular(test_dim2D_1, :Y, -0.5)
    @test isequal(parent(shift_2d_2), [2.5 5.5 8.5 NaN; 3.5 6.5 9.5 NaN; 4.5 7.5 10.5 NaN])

    shift_3d_1 = shift(test_dim3D_1, :Z, -0.5)
    @test isequal(
        parent(shift_3d_1[Y=At(5.0)]),
        [7.0 19.0 31.0 43.0 NaN; 8.0 20.0 32.0 44.0 NaN; 9.0 21.0 33.0 45.0 NaN],
    )
    shift_3d_2 = shift(test_dim3D_1, :Z, :Y, [-0.5, -1.0, -1.5, -2.0])
    @test all(isnan, shift_3d_2[Z=At(14)])
    shift_amount_along_y = DimArray([-0.5, -1.0, -1.5, -2.0], Y([5.0, 6.0, 7.0, 8]))
    @test isequal(parent(shift_3d_2), shift(test_dim3D_1, :Z, shift_amount_along_y))

    @test_throws ArgumentError shift(test_dim3D_1, :Z, test_dim2D_1)
    @test_throws ArgumentError shift(test_dim3D_1, :Z, :X, [-0.5, -1.0])
    @test_throws ArgumentError shift(test_dim3D_1, :W, :Y, [-0.5, -1.0, -1.5, -2.0])
    @test_throws ArgumentError shift(test_dim3D_1, :Z, :W, [-0.5, -1.0, -1.5, -2.0])
    @test_throws ArgumentError shift(shift_amount_along_y, :Y, :W, [-0.5])
end

function test_DimArray3D_irregular()
    x = range(1.0, stop = 3.0, length = 3)
    y = collect(range(5, stop = 8, length = 4))
    z = [10.0, 10.9, 12.0, 13.0, 14.0]
    data = 1.0:60.0 |> collect |> d->reshape(d, 3, 4, 5)
    da = DimArray(data, (X = x, Y = y, Z = z))
    return da
end

function test_DimArray3D_irregular_unordered()
    x = range(1.0, stop = 3.0, length = 3)
    y = collect(range(5, stop = 8, length = 4))
    z = [10.0, 12.0, 10.9, 13.0, 14.0]
    data = 1.0:60.0 |> collect |> d->reshape(d, 3, 4, 5)
    da = DimArray(data, (X = x, Y = y, Z = z))
    return da
end


@testset "Test for shift (2): Irregura DimArray" begin
    test_dim3D_2 = test_DimArray3D_irregular()

    shift_3d_2_2 = shift(test_dim3D_2, :Z, :Y, [-0.5, -1.0, -1.5, -2.0])
    @test all(isnan, shift_3d_2_2[Z=At(14)])
    shift_3d_2_3 = shift(test_dim3D_2, :Z, -0.5)

    test_dim3D_3 = test_DimArray3D_irregular_unordered()
    @test_throws ArgumentError shift(test_dim3D_3, :Z, :Y, [-0.5, -1.0, -1.5, -2.0])
    @test_throws ArgumentError shift(test_dim3D_3, :Z, -0.5)
end


@testset "Test for cat_arpes" begin
    data1 = test_ARPESData()
    data2 = test_ARPESData()
    data3 = test_ARPESData()
    data4 = test_spd_standard()

    data_1_1 = add_dim(data1, :w, 0.0)
    data_2_1 = add_dim(data2, :w, 1.0)
    data_3_1 = add_dim(data3, :w, 3.0)
    #    stack_1_to_3 = sort(cat([data_2_1, data_1_1, data_3_1]; dims = :w); dims = :w)
    #    @test size(stack_1_to_3) == (10, 10, 10, 3)
    #
    #    @test hasdim(stack_1_to_3, :w)
    #    @test lookup(stack_1_to_3, :w) == [0.0, 1.0, 3.0]
    #
    stack_1_to_3_alt = cat_arpes(data_2_1, data_1_1, data_3_1; dims = :w)
    #cat_arpes(data_2_1, data_1_1, data_3_1; dims = Dim{:w}([1.0, 0.0, 3.0]))
    @test parent(lookup(stack_1_to_3_alt, :w)) == [1.0, 0.0, 3.0]
    @test_throws ArgumentError cat_arpes(data_2_1, data_1_1, data3; dims = :w)
end

@testset "Test for rebin" begin
    spd_standard = test_spd_standard()
    # Rebin along phi with new edges
    rebin_spd = rebin(spd_standard, :phi, 10)
    rebin_spd_by_bins = rebin(spd_standard, :phi, Bins(10))
    @test size(rebin_spd) == (10, 601)
    @test parent(rebin_spd) == parent(rebin(spd_standard, :phi, Bins(10)))
    @test ARPES._is_equal_spacing(lookup(rebin_spd, :phi))
    @test metadata(rebin_spd)[:history][end] == "rebin by 10 along phi"
    @test metadata(rebin_spd)[:history][end] == "rebin by Bins(identity, 10) along phi"
end
