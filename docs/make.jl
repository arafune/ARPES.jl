using JuliARPES
using Documenter

DocMeta.setdocmeta!(JuliARPES, :DocTestSetup, :(using JuliARPES); recursive=true)

makedocs(;
    modules=[JuliARPES],
    authors="Ryuichi ARafune",
    sitename="JuliARPES.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
