using Test
using ARPES
using ARPES.IO: GenericLoader, SPDLoader, select_loader
using DimensionalData
using DimensionalData: name

@testset "load generic ITX fallback" begin
    path = tempname() * ".itx"
    try
        open(path, "w") do io
            write(
                io,
                """
                IGOR
                WAVES/N=4 testwave
                BEGIN
                1 2 3 4
                END
                X SetScale/P x, 0, 1, "", "testwave"
                """,
            )
        end
        @test select_loader(path, nothing) == GenericLoader

        data = load(path)
        @test data isa ARPESData
        @test name(dims(data)[1]) == :x
    finally
        rm(path; force = true)
    end
end

@testset "load errors are explicit" begin
    path = tempname() * ".foo"
    try
        write(path, "dummy")
        @test_throws ArgumentError load(path)
    finally
        rm(path; force = true)
    end

    path = tempname() * ".sp2"
    try
        write(path, "")
        @test select_loader(path, nothing) == SPDLoader
        @test_throws ArgumentError load(path; loc = "SPD")
    finally
        rm(path; force = true)
    end
end
