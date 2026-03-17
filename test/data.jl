using Test
using ARPES
using ARPES: add_dim, cat_arpes
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using ARPES: _is_equal_spacing
using DimensionalData: Dim, hasdim, lookup

if !@isdefined(test_ARPESData)
    include("fixture.jl")
end

@testset "test for shift (1)" begin
    test_dim2D_1 = test_DimArray2D()
    test_dim3D_1 = test_DimArray3D()
    shift_2d_1 = shift(test_dim2D_1, :Y, -0.5)
    @test isequal(parent(shift_2d_1), [2.5 5.5 8.5 NaN; 3.5 6.5 9.5 NaN; 4.5 7.5 10.5 NaN])

end


@testset "Test for vstack_arpes" begin
    # Originally, stack_arpes function was designed to stack ARPESData objects
    # along a new dimension, but the functionality I expected can be easily realized
    # by using the combination cat with sort. This test demonstrates how to do it.

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

@testset "Test for _is_equal_spacing" begin
    @test _is_equal_spacing([1])
    @test _is_equal_spacing([1, 2, 3])
    @test _is_equal_spacing([1.0, 2.0, 3.0, 4.0, 5.0])
    @test !_is_equal_spacing([1, 2, 4])
    @test !_is_equal_spacing([1.0, 2.0, 4.0])
    @test _is_equal_spacing([1, 2]; atol = 0.5)
end

