using MeasureBase
using Documenter

DocMeta.setdocmeta!(MeasureBase, :DocTestSetup, :(using MeasureBase); recursive=true)

makedocs(;
    modules=[MeasureBase],
    authors="Chad Scherrer <chad.scherrer@gmail.com> and contributors",
    repo="https://github.com/cscherrer/MeasureBase.jl/blob/{commit}{path}#{line}",
    sitename="MeasureBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cscherrer.github.io/MeasureBase.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cscherrer/MeasureBase.jl",
)
