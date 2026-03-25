using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using ARPES: _step
using DimensionalData: hasdim, name, lookup

if !@isdefined(test_ARPESData)
    include("fixture.jl")
end

@testset "Test for convert_dim" begin
    common_metadata = Dict(:unit => "nm", :description => "Test dimension")
    dim_range = Dim{:x}(0:0.5:10, metadata = common_metadata)
    dim_vector = Dim{:y}(collect(0:0.1:10), metadata = common_metadata)
    dim_vector_non_uniform =
        Dim{:z}([0, 0.1, 0.3, 0.6, 1.0], metadata = Dict(:unit => "nm"))

    linear_func = x -> 2x + 1.0
    nonlinear_func = x -> x^2 + 3x + 2.0
    @testset "Test by range based Dimension" begin
        @testset "Test for convert_dim by linear function" begin
            converted_range = convert_dim(dim_range, linear_func)
            @test converted_range isa DimensionalData.Dimension
            @test parent(lookup(converted_range)) isa AbstractRange
            @test step(converted_range) == 2 * step(dim_range)
            @test metadata(converted_range) == common_metadata
        end
        @testset "Test for convert_dim by non-linear function" begin
            converted_range_nonlinear = convert_dim(dim_range, nonlinear_func)
            @test converted_range_nonlinear isa DimensionalData.Dimension
            @test parent(lookup(converted_range_nonlinear)) isa AbstractVector
            @test metadata(converted_range_nonlinear) == common_metadata
            @test length(converted_range_nonlinear) == length(dim_range)
            @test_throws MethodError step(converted_range_nonlinear)
        end
    end
    @testset "Test by vector based Dimension" begin
        @testset "Test for convert_dim by linear function" begin
            converted_vector = convert_dim(dim_vector, linear_func)
            @test converted_vector isa DimensionalData.Dimension
            @test parent(lookup(converted_vector)) isa AbstractRange
            @test metadata(converted_vector) == common_metadata
            @test length(converted_vector) == length(dim_vector)
            @test step(converted_vector) == 2 * _step(dim_vector)
        end
        @testset "Test for convert_dim by non-linear function" begin
            converted_vector = convert_dim(dim_vector, nonlinear_func)
            @test converted_vector isa DimensionalData.Dimension
            @test parent(lookup(converted_vector)) isa AbstractVector
            @test metadata(converted_vector) == common_metadata
            @test length(converted_vector) == length(dim_vector)
            @test_throws MethodError step(converted_vector)
        end
    end
    @testset "Test by non-uniform vector based Dimension" begin
        @testset "Test for convert_dim by linear function" begin
            converted_non_uniform = convert_dim(dim_vector_non_uniform, linear_func)
            @test converted_non_uniform isa DimensionalData.Dimension
            @test parent(lookup(converted_non_uniform)) isa AbstractVector
            @test metadata(converted_non_uniform) == Dict(:unit => "nm")
            @test length(converted_non_uniform) == length(dim_vector_non_uniform)
            @test_throws MethodError step(converted_non_uniform)
        end
        @testset "Test for convert_dim by non-linear function" begin
            converted_non_uniform = convert_dim(dim_vector_non_uniform, nonlinear_func)
            @test converted_non_uniform isa DimensionalData.Dimension
            @test parent(lookup(converted_non_uniform)) isa AbstractVector
            @test metadata(converted_non_uniform) == Dict(:unit => "nm")
            @test length(converted_non_uniform) == length(dim_vector_non_uniform)
            @test_throws MethodError step(converted_non_uniform)
        end
    end
    @testset "Test for error handling in convert_dim" begin
        @test_throws ArgumentError convert_dim(dim_range, 1.0)
        @test_throws ArgumentError convert_dim(dim_vector, "not a function")
        @test_throws ArgumentError convert_dim(dim_vector_non_uniform, missing)
    end
end

@testset "Test for convert_dim method dispatch" begin
    common_metadata = Dict(:unit => "nm", :description => "Test dimension")
    dim_range = Dim{:x}(0:0.5:10, metadata = common_metadata)
    linear_func = x -> 2x + 1.0

    @testset "Default label overload" begin
        converted = convert_dim(dim_range, linear_func)
        @test converted isa DimensionalData.Dimension
        @test name(converted) == :x
        @test parent(lookup(converted)) isa AbstractRange
        @test metadata(converted) == common_metadata
        @test step(converted) == 2 * step(dim_range)
    end

    @testset "Non-function overload" begin
        @test_throws ArgumentError convert_dim(dim_range, 1.0)
    end
end

@testset "Test for convert_dim on AbstractDimArray" begin
    data = test_ARPESData()
    linear_func = x -> 2x + 1.0
    nonlinear_func = x -> x^2 + 3x + 2.0

    @testset "Symbol selector" begin
        converted = convert_dim(data, :phi, linear_func)
        @test converted isa ARPESData
        @test name.(dims(converted)) == name.(dims(data))
        @test name(dims(converted)[1]) == :phi
        @test parent(lookup(dims(converted)[1])) isa AbstractRange
        @test step(dims(converted)[1]) == 2 * step(dims(data)[1])
        @test dims(converted)[2] == dims(data)[2]
        @test metadata(converted) == metadata(data)
    end

    @testset "Dimension selector" begin
        converted = convert_dim(data, dims(data)[2], nonlinear_func)
        @test converted isa ARPESData
        @test name.(dims(converted)) == name.(dims(data))
        @test name(dims(converted)[2]) == name(dims(data)[2])
        @test parent(lookup(dims(converted)[2])) isa AbstractVector
        @test_throws MethodError step(dims(converted)[2])
        @test dims(converted)[1] == dims(data)[1]
        @test metadata(converted) == metadata(data)
    end

    @testset "Error handling" begin
        @test_throws ArgumentError convert_dim(data, :delay, linear_func)
        @test_throws ArgumentError convert_dim(data, :phi, 1.0)
    end
end

@testset "Test for shift_dim" begin
    x = range(1, stop = 3, length = 10)
    y = collect(range(-1, stop = 8, length = 10))
    z = 1:10
    data = rand(10, 10, 10)
    da = DimArray(data, (X(x), Y(y), Z(z)))
    x_dim = da.dims[1]
    y_dim = da.dims[2]

    dim_x_shifted = shift_dim(x_dim, 5)
    @test dim_x_shifted isa DimensionalData.Dimension
    @test dim_x_shifted[1] == 6
    @test dim_x_shifted[end] == 8
    @test length(dim_x_shifted) == length(x_dim)

    dim_y_shifted = shift_dim(y_dim, -2.2)
    @test dim_y_shifted isa DimensionalData.Dimension
    @test dim_y_shifted[1] == -3.2
    @test dim_y_shifted[end] == 5.8
    @test length(dim_y_shifted) == length(y_dim)
end

@testset "Test for negate_dim" begin
    x = range(1, stop = 3, length = 10)
    y = collect(range(-1, stop = 8, length = 10))
    z = 1:10
    data = rand(10, 10, 10)
    da = DimArray(data, (X(x), Y(y), Z(z)))
    x_dim = da.dims[1]
    y_dim = da.dims[2]

    dim_x_negated = negate_dim(x_dim)
    @test dim_x_negated isa DimensionalData.Dimension
    @test dim_x_negated[1] == -1
    @test dim_x_negated[end] == -3
    @test length(dim_x_negated) == length(x_dim)

    dim_y_negated = negate_dim(y_dim)
    @test dim_y_negated isa DimensionalData.Dimension
    @test dim_y_negated[1] == 1
    @test dim_y_negated[end] == -8
    @test length(dim_y_negated) == length(y_dim)
end

@testset "Test for shift_dim on ARPESData" begin
    data = test_ARPESData()
    shifted_phi_data = shift_dim(data, :phi, 5.0)
    @test shifted_phi_data isa ARPESData
    #    @test shifted_phi_data.dims[1][1] == 6
    #    @test shifted_phi_data.dims[1][end] == 15
    #    @test shifted_phi_data.dims[2] == data.dims[2]
    #
    shifted_phi_data2 = shift_dim(data, dims(data)[1], 5.0)
    @test shifted_phi_data2 isa ARPESData
    shifted_eV_data = shift_dim(data, eV, -0.5)
    @test shifted_eV_data isa ARPESData
    shifted_both = shift_dim(data, :eV, -0.5, phi, 5)
    @test shifted_both isa ARPESData
    @test_throws ArgumentError shift_dim(data, :eV, 0.1, :ky)
    @test_throws ArgumentError shift_dim(data, :delay, 0.1)
    @test_throws ArgumentError shift_dim(data, ch2, 0.1)
end


@testset "Test for add_dim" begin
    data = test_ARPESData()
    data_4d = add_dim(data, :w, 0.0)
    @test size(data_4d) == (10, 10, 10, 1)
    @test hasdim(data_4d, :w) == true
    data_5d = add_dim(data_4d, :β)
    @test size(data_5d) == (10, 10, 10, 1, 1)
    @test hasdim(data_5d, :β) == true
    data_6d = add_dim(data_5d, "δ", 0.0)
    @test size(data_6d) == (10, 10, 10, 1, 1, 1)
    @test hasdim(data_6d, :δ) == true
    @test_throws ArgumentError add_dim(data, :phi, 1.0)
    @test_throws ArgumentError add_dim(data, :χ)
end

@testset "Test for rename_dim" begin
    data = test_ARPESData()
    renamed_data = rename_dim(data, :phi => :angle)
    @test hasdim(renamed_data, :angle) == true
    @test hasdim(renamed_data, :phi) == false
    @test_throws ArgumentError rename_dim(data, :delay => :time)
end
