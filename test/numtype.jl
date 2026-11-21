# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: StdNormal, StdUniform, LebesgueBase, Lebesgue, Dirac
using MeasureBase: logdensity_def, logdensity_rel, weightedmeasure, superpose, logdensities, mbind
using Static: static
import ForwardDiff

@testset "number type of log-densities" begin
    x = 0.5f0
    xf = randn(Float32, 3)
    uf = rand(Float32, 3)
    mix = superpose(weightedmeasure(log(0.3f0), StdNormal()), weightedmeasure(log(0.7f0), StdUniform()))

    @test @inferred(logdensityof(StdNormal(), x)) isa Float32
    @test @inferred(logdensityof(StdNormal()^3, xf)) isa Float32
    @test @inferred(MeasureBase.unsafe_logdensityof(StdNormal()^3, xf)) isa Float32
    @test @inferred(logdensity_def(LebesgueBase()^3, xf)) isa Float32
    @test @inferred(logdensity_def(LebesgueBase(), x)) isa Float32
    @test @inferred(logdensity_rel(StdNormal(), StdUniform(), x)) isa Float32
    @test @inferred(logdensity_rel(StdNormal()^3, StdUniform()^3, uf)) isa Float32
    @test @inferred(logdensity_rel(StdNormal()^3, StdUniform()^3, xf)) isa Float32
    @test @inferred(logdensity_rel(Lebesgue(), Dirac(1f0), 2f0)) isa Float32
    @test @inferred(logdensity_rel(Dirac(1f0), Lebesgue(), 1f0)) isa Float32
    @test @inferred(logdensityof(Dirac(1f0), 1f0)) isa Float32
    @test @inferred(logdensityof(weightedmeasure(0.5f0, StdNormal())^3, xf)) isa Float32
    @test @inferred(logdensityof(weightedmeasure(static(0.5), StdNormal())^3, xf)) isa Float32
    @test @inferred(logdensityof(weightedmeasure(0.5, StdNormal()), x)) isa Float32
    @test @inferred(logdensityof(mix, x)) isa Float32
    @test @inferred(logdensityof(mix^3, uf)) isa Float32
    @test @inferred(logdensities(StdNormal(), randn(Float32, 4))) isa Vector{Float32}
    @test @inferred(logdensities(StdNormal()^3, randn(Float32, 3, 4))) isa Vector{Float32}

    # Weights that carry derivatives keep them:
    @test ForwardDiff.derivative(w -> logdensityof(weightedmeasure(w, StdNormal()), 0.3), 0.1) ≈ 1
    @test ForwardDiff.derivative(w -> logdensityof(weightedmeasure(w, StdNormal())^2, [0.3, 0.1]), 0.1) ≈ 2
end
