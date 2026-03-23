using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using ARPES: _apply_along_dim
using DimensionalData
using DimensionalData: Dim, hasdim, lookup
using Random
using Statistics


@testset "_apply_along_dim" begin

    # -------------------------------
    # 1. 1D: identity (no-op) test
    # Ensures data and metadata are preserved
    # -------------------------------
    t = range(0, 10, length = 50)
    A = DimArray(collect(t), (Dim{:t}(t),))

    B = _apply_along_dim(A, :t, x -> x)

    @test parent(B) == parent(A)
    @test dims(B) == dims(A)
    @test size(B) == size(A)

    # -------------------------------
    # 2. 1D: simple transformation
    # Verify function is correctly applied
    # -------------------------------
    B = _apply_along_dim(A, :t, x -> x .+ 1)

    @test parent(B) ≈ parent(A) .+ 1

    # -------------------------------
    # 3. 2D: apply along first dimension
    # Each column should be processed independently
    # -------------------------------
    x = 1:4
    y = 1:3
    data = reshape(collect(1.0:12.0), 4, 3)

    A2 = DimArray(data, (Dim{:x}(x), Dim{:y}(y)))

    # Subtract mean along :x (column-wise)
    B2 = _apply_along_dim(A2, :x, v -> v .- mean(v))

    for j in axes(data, 2)
        @test mean(parent(B2)[:, j]) ≈ 0
    end

    @test size(B2) == size(A2)
    @test dims(B2) == dims(A2)

    # -------------------------------
    # 4. Apply along second dimension
    # Each row should be processed independently
    # -------------------------------
    B3 = _apply_along_dim(A2, :y, v -> v .- mean(v))

    for i in axes(data, 1)
        @test mean(parent(B3)[i, :]) ≈ 0
    end

    # -------------------------------
    # 5. Position-dependent function
    # Ensures correct indexing inside each slice
    # -------------------------------
    B4 = _apply_along_dim(A2, :x, v -> v .* (1:length(v)))

    for j in axes(data, 2)
        expected = data[:, j] .* (1:length(x))
        @test parent(B4)[:, j] == expected
    end

    # -------------------------------
    # 6. Type stability
    # Output element type should be preserved
    # -------------------------------
    Afloat = DimArray(rand(Float32, 5, 5), (Dim{:x}(1:5), Dim{:y}(1:5)))
    Bfloat = _apply_along_dim(Afloat, :x, x -> x .* 2)

    @test eltype(Bfloat) == Float32

    # -------------------------------
    # 7. NaN propagation
    # Function should not silently remove NaNs
    # -------------------------------
    A_nan = DimArray([1.0, NaN, 3.0], (Dim{:t}(1:3),))
    B_nan = _apply_along_dim(A_nan, :t, x -> x)

    @test isnan(parent(B_nan)[2])

    # -------------------------------
    # 8. 3D test
    # Ensure correct behavior in higher dimensions
    # -------------------------------
    A3 = DimArray(rand(4, 3, 2), (Dim{:x}(1:4), Dim{:y}(1:3), Dim{:z}(1:2)))

    B3d = _apply_along_dim(A3, :x, v -> v .- mean(v))

    for j = 1:3, k = 1:2
        @test mean(parent(B3d)[:, j, k]) ≈ 0 atol=1e-12
    end

    @test size(B3d) == size(A3)
    @test dims(B3d) == dims(A3)

    # -------------------------------
    # 9. Consistency with mapslices
    # Compare with Base.mapslices behavior
    # -------------------------------
    f(v) = v .- mean(v)

    B_ref = mapslices(f, parent(A2); dims = 1)  # dims=1 corresponds to :x
    B_test = _apply_along_dim(A2, :x, f)

    @test parent(B_test) ≈ B_ref

end
