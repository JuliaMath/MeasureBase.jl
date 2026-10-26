# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

import Distributions
import MeasureBase
import AffineMaps
using LinearAlgebra: I

@testset "Measure interface" begin
    c0 = Distributions.Weibull(0.7, 1.3)
    c1 = Distributions.MvNormal([0.7, 0.9], [1.4 0.5; 0.5 1.1])

    d0 = Distributions.Poisson(0.7)
    d1 = Distributions.product_distribution(Distributions.Poisson.([0.7, 1.4]))

    for d in [c0, c1, d0, d1]
        x = rand(d)
        @test @inferred(MeasureBase.logdensity_def(d, x)) == Distributions.logpdf(d, x)
        @test @inferred(MeasureBase.unsafe_logdensityof(d, x)) == Distributions.logpdf(d, x)

        MeasureBase.Interface.test_interface(d)
    end

    @test @inferred(MeasureBase.basemeasure(c0)) == MeasureBase.Lebesgue(MeasureBase.ℝ)
    @test @inferred(MeasureBase.basemeasure(c1)) == MeasureBase.Lebesgue(MeasureBase.ℝ) ^ 2

    @test @inferred(MeasureBase.insupportof(c0, 3)) == true
    @test @inferred(MeasureBase.insupportof(c0, -3)) == false
    @test @inferred(MeasureBase.insupportof(c1, [0.1, 0.2])) == true
    @test @inferred(MeasureBase.insupportof(d0, 3)) == true
    @test @inferred(MeasureBase.insupportof(d0, 3.2)) == false
    @test @inferred(MeasureBase.insupportof(d1, [1, 2])) == true
    @test @inferred(MeasureBase.insupportof(d1, [1.1, 2.2])) == false

    @test MeasureBase.paramnames(c0) == (:α, :θ)
    if VERSION >= v"1.8"
        @test @inferred(MeasureBase.paramsof(c0)) == (α = 0.7, θ = 1.3)
    else
        # v1.6 can't type-infer this:
        @test (MeasureBase.paramsof(c0)) == (α = 0.7, θ = 1.3)
    end

    @test MeasureBase.mintegrate(x -> Distributions.Normal(x, 0), Distributions.Normal()) isa MeasureBase.DensityMeasure
end


# Needs AffineMaps >= 0.3.4 (earlier versions give per-column Jacobians of
# scalar factors and scalar ones for Add):
@testset "affine pushforward densities" begin
    b = [1.0, 2.0, 3.0]
    x = [0.3, -1.2, 2.0]
    X = [0.3 1.0; -1.2 0.5; 2.0 -0.7]
    for (f, d) in [
        (AffineMaps.Mul(2.0), Distributions.MvNormal(zeros(3), 4.0 * I(3))),
        (AffineMaps.Add(b), Distributions.MvNormal(b, 1.0 * I(3))),
        (AffineMaps.MulAdd(2.0, b), Distributions.MvNormal(b, 4.0 * I(3))),
    ]
        ν = MeasureBase.pushfwd(f, MeasureBase.StdNormal()^3)
        @test MeasureBase.logdensityof(ν, x) ≈ Distributions.logpdf(d, x)
        @test MeasureBase.logdensities(ν, X) ≈ Distributions.logpdf(d, X)
    end
end
