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


@testset "bare distributions behave like their measures" begin
    for d in (Distributions.Uniform(1.0, 2.0), Distributions.Dirichlet([1.0, 2.0, 3.0]), Distributions.Categorical([0.2, 0.5, 0.3]),
              Distributions.Multinomial(3, [0.2, 0.8]), Distributions.Poisson(2.0), Distributions.product_distribution(Distributions.Uniform.([1.0, 2.0], [2.0, 3.0])))
        @test Distributions.insupport(d, MeasureBase.testvalue(d))
        @test MeasureBase.basemeasure(d) == MeasureBase.basemeasure(MeasureBase.asmeasure(d))
    end
    @test MeasureBase.basemeasure(Distributions.Dirichlet([1.0, 2.0, 3.0])) == MeasureBase.Lebesgue(MeasureBase.StdSimplex(3))
    @test MeasureBase.basemeasure(Distributions.Poisson(2.0)) == MeasureBase.Counting(MeasureBase.ℤ)
    @test !MeasureBase.insupportof(Distributions.Normal(), Inf)
    # Constructors that took bare distributions:
    @test MeasureBase.logdensityof(MeasureBase.weightedmeasure(0.5, Distributions.Normal()), 0.3) ≈ Distributions.logpdf(Distributions.Normal(), 0.3) + 0.5
    p = MeasureBase.pushfwd(exp, Distributions.Uniform(1.0, 2.0))
    @test exp(1.0) <= MeasureBase.testvalue(p) <= exp(2.0)
    @test MeasureBase.logdensityof(p, exp(1.5)) ≈ -1.5
end

@testset "parameters of distributions" begin
    I2 = Matrix(1.0I, 2, 2)
    for (d, names) in [
        (Distributions.Weibull(0.7, 1.3), (:α, :θ)), (Distributions.Dirichlet([1.0, 2.0, 3.0]), (:alpha,)),
        (Distributions.Categorical([0.2, 0.8]), (:p,)), (Distributions.DiscreteUniform(1, 6), (:a, :b)),
        (Distributions.DiscreteNonParametric([0.5, 1.5], [0.4, 0.6]), (:support, :p)),
        (Distributions.Wishart(3.0, I2), (:df, :S)), (Distributions.LKJ(2, 1.0), (:d, :η)),
        (Distributions.MvTDist(3.0, [0.0, 1.0], I2), (:df, :μ, :Σ)), (Distributions.MvNormal([0.0, 1.0], I2), (:μ, :Σ)),
        (Distributions.truncated(Distributions.Normal(), 0.0, 1.0), (:p1, :p2, :p3, :p4)),
        (Distributions.censored(Distributions.Normal(), 0.0, 1.0), (:p1, :p2, :p3, :p4)),
    ]
        m = MeasureBase.asmeasure(d)
        @test MeasureBase.paramnames(d) == MeasureBase.paramnames(m) == names
        @test values(MeasureBase.paramsof(d)) == values(MeasureBase.paramsof(m)) == Distributions.params(d)
    end
    @test @inferred(MeasureBase.paramsof(Distributions.Dirichlet([1.0, 2.0]))) == (alpha = [1.0, 2.0],)
    @test MeasureBase.paramsof(MeasureBase.asmeasure(Distributions.Dirichlet([1.0, 2.0]))^3) == (alpha = [1.0, 2.0],)
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
