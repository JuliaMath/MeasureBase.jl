# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using MeasureBase

# Doctest setup
DocMeta.setdocmeta!(
    MeasureBase,
    :DocTestSetup,
    :(using MeasureBase);
    recursive=true,
)

makedocs(
    sitename = "MeasureBase",
    modules = [MeasureBase],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://juliamath.github.io/MeasureBase.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/juliamath/MeasureBase.jl.git",
    forcepush = true,
    push_preview = true,
)
