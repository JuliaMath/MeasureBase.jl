# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using MeasureBase

# Doctest setup
DocMeta.setdocmeta!(MeasureBase, :DocTestSetup, :(using MeasureBase); recursive = true)

# Sorts docstrings into the sections of api.md. Functions of Base and other
# packages that MeasureBase extends count as exported API:
module MeasureBaseDocs
    using MeasureBase
    function api_section(obj)
        obj isa Union{Function,DataType,UnionAll,Module} || return :internal
        name = nameof(obj)
        if Base.isexported(MeasureBase, name) || !(parentmodule(obj) in (MeasureBase, MeasureBase.MeasureOperators))
            :exported
        elseif Base.ispublic(MeasureBase, name)
            :public
        else
            :internal
        end
    end
end

makedocs(
    sitename = "MeasureBase",
    modules = [MeasureBase, MeasureBase.MeasureOperators],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://juliamath.github.io/MeasureBase.jl/stable/",
    ),
    pages = [
        "Home" => "index.md",
        "Concepts" => "concepts.md",
        "API" => "api.md",
        "Internals" => "internals.md",
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
