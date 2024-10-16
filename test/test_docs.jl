# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test
using MeasureBase
import Documenter

Documenter.DocMeta.setdocmeta!(
    MeasureBase,
    :DocTestSetup,
    :(using MeasureBase);
    recursive = true,
)
Documenter.doctest(MeasureBase)
