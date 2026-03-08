using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using DimensionalData

@testset "test energy_definition in eV metadata" begin
    @test metadata(dims(spd_standard)[dimnum(spd_standard, :eV)])[:energy_definition] ==
          FinalStateEnergy
end



@testset "Test for shift_dim" begin
    x = range(1, stop = 3, length = 10)
    y = 1:10
    z = 1:10
    data = rand(10, 10, 10)
    da = DimArray(data, (X(x), Y(y), Z(z)))
    x_dim = da.dims[1]

    dim_x_shifted = shift_dim(x_dim, 5)
    @test dim_x_shifted isa DimensionalData.Dimension
    @test dim_x_shifted[1] == 6
    @test dim_x_shifted[end] == 8
    @test length(dim_x_shifted) == length(x_dim)
end


