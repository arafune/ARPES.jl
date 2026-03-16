using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using DimensionalData: Dim, hasdim

if !@isdefined(test_ARPESData)
    include("fixture.jl")
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
    @test size(data_4d)== (10, 10, 10, 1)
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
