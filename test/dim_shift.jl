using Test
using DimensionalData
using ARPES


@testset "Test for shift_dim" begin
    x = range(1:3, length = 10)
    y = 1:10
    z = 1:10
    data = rand(10, 10, 10)
    da = DimArray(data, (X(x), Y(y), Z(z)))
    x_dim = da.dims[1]

    dim_x_shifted = shift_dim(x_dim, 5)
    @test dim_x_shifted isa DimensionalData.Dimension
    @test dim_x_shifted[1] == 6
    @teest dim_x_shifted[end] == 8
    @test length(dim_x_shifted) == length(x_dim)
end


