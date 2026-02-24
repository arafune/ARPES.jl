using ARPES
using Documenter

DocMeta.setdocmeta!(ARPES, :DocTestSetup, :(using ARPES); recursive = true)

makedocs(;
    modules = [ARPES],
    authors = "Ryuichi Arafune",
    sitename = "ARPES.jl",
    format = Documenter.HTML(; edit_link = "main", assets = String[]),
    pages = ["Home" => "index.md"],
)
