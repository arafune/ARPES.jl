using Test
using Interpolations

@testset "fast_interpolate tests" begin

    # --- 1. 2D Regular Grid (BSpline) ---
    @testset "2D Regular" begin
        c0, c1 = 0.0:1.0:2.0, 0.0:1.0:2.0
        values = [
            1.0 2.0 3.0;
            4.0 5.0 6.0;
            7.0 8.0 9.0
        ]

        # Test point exactly on a node
        @test fast_interpolate([0.0], [0.0], c0, c1, values) ≈ [1.0]
        # Test point in the middle
        @test fast_interpolate([0.5], [0.5], c0, c1, values) ≈ [3.0]
        # Test out of bounds (NaN)
        @test isnan(fast_interpolate([-1.0], [0.0], c0, c1, values)[1])
    end

    # --- 2. 2D Irregular Grid (Gridded) ---
    @testset "2D Irregular" begin
        c0, c1 = [0.0, 1.0, 5.0], [0.0, 2.0, 10.0] # Non-uniform
        values = ones(3, 3)
        p0, p1 = [0.5 5.0; -1.0 2.5], [1.0 10.0; 0.0 5.0] # 2x2 matrix input

        res = fast_interpolate(p0, p1, c0, c1, values)
        @test size(res) == (2, 2)
        @test res[1, 1] ≈ 1.0
        @test isnan(res[2, 1]) # p0 = -1.0 is out
    end

    # --- 3. 3D Regular Grid (BSpline) ---
    @testset "3D Regular" begin
        c0 = c1 = c2 = 0.0:1.0:1.0
        values = ones(2, 2, 2)
        @test fast_interpolate([0.5], [0.5], [0.5], c0, c1, c2, values) ≈ [1.0]
    end

    # --- 4. 3D Irregular Grid (Gridded) ---
    @testset "3D Irregular" begin
        c0, c1, c2 = [0.0, 1.0, 10.0], [0.0, 1.0], [0.0, 1.0]
        values = zeros(3, 2, 2)
        values[2, 1, 1] = 10.0 # point at (1.0, 0.0, 0.0)

        @test fast_interpolate([1.0], [0.0], [0.0], c0, c1, c2, values) ≈ [10.0]
        @test isnan(fast_interpolate([11.0], [0.0], [0.0], c0, c1, c2, values)[1])
    end
end
