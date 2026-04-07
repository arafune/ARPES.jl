using Test
using ARPES
using DimensionalData

@testset "read_itx basic test" begin
    file_path = joinpath(pkgdir(ARPES), "testdata", "arpes_ch_resolved.itx")
    data = ARPES.IO.Format.read_itx(file_path)
    meta_data = metadata(data)
    @test haskey(meta_data, :created_date)
    @test haskey(meta_data, :created_by)
    @test haskey(meta_data, :scan_mode)
    @test haskey(meta_data, :number_of_scans)
    @test haskey(meta_data, :excitation_energy)
    @test haskey(meta_data, :user_comment)
    @test haskey(meta_data, :x_scale)
end
