using ARPES
using Test
using Coverage

include("./fixture.jl")

@testset "ARPES.jl" begin
    @testset "function for Dict as metadata" begin
        include("dict.jl")
    end
    @testset "Core" begin
        include("utils.jl")
        include("dims.jl")
        include("transform.jl")
        include("filter.jl")
        include("derivative.jl")
        include("stitch.jl")
    end
    @testset "k_conv" begin
        include("k_conv/preprocess.jl")
        include("k_conv/interpolation.jl")
        include("k_conv/k_conversion.jl")
    end
    @testset "IO" begin
        @testset "Location" begin
            include("io/location/location.jl")
        end
        @testset "FORMATS" begin
            include("io/formats/itx.jl")
        end
        @testset "registry" begin
            include("io/registry.jl")
        end
    end
end

testdata_dir = joinpath(@__DIR__, "..", "testdata")

if isdir(testdata_dir)
    include(joinpath(@__DIR__, "fixture_testdata.jl"))
    @testset "ARPESData with test data" begin
        include("ARPESData_testdata.jl")
    end
    @testset "Core" begin
        include("transform.jl")
    end
    @testset "k_conv" begin
        include("k_conv/preprocess_testdata.jl")
        include("k_conv/k_conversion_testdata.jl")
    end
    @testset "IO basic" begin
        include("io/io_testdata.jl")
    end
    @testset "Location" begin
        include("io/location/spd_testdata.jl")
    end
    @testset "FORMATS" begin
        include("io/formats/itx_testdata.jl")
    end

    @testset "registry" begin
        include("io/registry_testdata.jl")
    end
end
