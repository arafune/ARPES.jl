using ARPES
using Coverage
using Test

file_path1 = joinpath(pkgdir(ARPES), "testdata", "arpes_ch_resolved.itx")
file_path2 = joinpath(pkgdir(ARPES), "testdata", "spd_standard.itx")

meta_data = Dict(:beta=>0.0)
arpes_ch_resolved = load(file_path1, loc = "SPD", extra_metadata = meta_data)  ## Never change in the test
spd_standard = load(file_path2, loc = "SPD")                                   ## Never change in the test

@testset "ARPES.jl" begin
    @testset "Core" begin
        include("ARPESData.jl")
        include("dims.jl")
    end
    @testset "IO" begin
        @testset "IO basic" begin
            include("io/io.jl")
        end
        @testset "Location" begin
            include("io/location/location.jl")
            include("io/location/spd.jl")
        end
        @testset "FORMATS" begin
            include("io/formats/itx.jl")
        end
        @testset "registry" begin
            include("io/registry.jl")
        end
    end
end
