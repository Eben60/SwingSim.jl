using SwingSim
using Documenter

DocMeta.setdocmeta!(SwingSim, :DocTestSetup, :(using SwingSim); recursive=true)

makedocs(;
    modules=[SwingSim],
    authors="Eben60",
    sitename="SwingSim.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
