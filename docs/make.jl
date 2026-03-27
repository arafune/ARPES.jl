using ARPES
using Documenter

DocMeta.setdocmeta!(ARPES, :DocTestSetup, :(using ARPES); recursive = true)

makedocs(;
    modules = [ARPES],
    authors = "Ryuichi Arafune",
    sitename = "ARPES.jl",
    format = Documenter.HTML(; edit_link = "main", assets = String[]),
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
        "User Guide" => [
            "Working with ARPESData" => "guide/data-model.md",
            "Data loading" => "guide/io.md",
            "Analysis utilities" => "guide/analysis.md",
            "k-space conversion" => "guide/k-space.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/arafune/ARPES.jl", push_preview = true)
