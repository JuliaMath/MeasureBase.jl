# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

import Distributions
using Distributions: Distribution
import MeasureBase
using MeasureBase: AbstractMeasure

@testset "Measure interface" begin
    d = Distributions.Weibull()
    @test @inferred(AbstractMeasure(d)) isa AbstractMeasure
    @test @inferred(AbstractMeasure(d)) isa DistributionMeasure
    @test @inferred(convert(AbstractMeasure, d)) isa AbstractMeasure
    @test @inferred(convert(AbstractMeasure, d)) isa DistributionMeasure
    @test @inferred(Distribution(AbstractMeasure(d))) === d
    @test @inferred(convert(Distribution, convert(AbstractMeasure, d))) === d


    c0 = AbstractMeasure(Distributions.Weibull(0.7, 1.3))
    c1 = AbstractMeasure(Distributions.MvNormal([0.7, 0.9], [1.4 0.5; 0.5 1.1]))

    d0 = AbstractMeasure(Distributions.Poisson(0.7))
    d1 = AbstractMeasure(Distributions.product_distribution(Distributions.Poisson.([0.7, 1.4])))

    for μ in [c0, c1, d0, d1]
        d = Distribution(μ)
        x = rand(μ)
        @test @inferred(MeasureBase.logdensity_def(μ, x)) == Distributions.logpdf(d, x)
        @test @inferred(MeasureBase.unsafe_logdensityof(μ, x)) == Distributions.logpdf(d, x)

        MeasureBase.Interface.test_interface(d)
    end

    @test @inferred(MeasureBase.basemeasure(c0)) == MeasureBase.Lebesgue(MeasureBase.ℝ)
    @test @inferred(MeasureBase.basemeasure(c1)) == MeasureBase.Lebesgue(MeasureBase.ℝ) ^ 2

    @test @inferred(MeasureBase.insupport(c0, 3)) == true
    @test @inferred(MeasureBase.insupport(c0, -3)) == false
    @test @inferred(MeasureBase.insupport(c1, [0.1, 0.2])) == true
    @test @inferred(MeasureBase.insupport(d0, 3)) == true
    @test @inferred(MeasureBase.insupport(d0, 3.2)) == false
    @test @inferred(MeasureBase.insupport(d1, [1, 2])) == true
    @test @inferred(MeasureBase.insupport(d1, [1.1, 2.2])) == false

    @test MeasureBase.paramnames(c0) == (:α, :θ)
    if VERSION >= v"1.8"
        @test @inferred(MeasureBase.params(c0)) == (α = 0.7, θ = 1.3)
    else
        # v1.6 can't type-infer this:
        @test (MeasureBase.params(c0)) == (α = 0.7, θ = 1.3)
    end
end

@testset "pushforwards of discrete distributions" begin
    using Distributions: Poisson, Multinomial, logpdf
    using AffineMaps: Mul
    m = AbstractMeasure(Poisson(2.0))
    p = MeasureBase.pushfwd(Mul(2.0), m)
    @test MeasureBase.logdensityof(p, 4.0) == MeasureBase.logdensityof(m, 2.0)
    @test MeasureBase.logdensityof(p, 3.0) == -Inf
    @test sum(y -> exp(MeasureBase.logdensityof(p, y)), 0.0:0.5:80.0) ≈ 1
    @test MeasureBase.logdensities(p, [0.0, 3.0, 4.0]) ≈ logpdf.(Poisson(2.0), [0.0, 1.5, 2.0])
    @test MeasureBase.logdensity_rel(p, AbstractMeasure(Poisson(4.0)), 4.0) ≈ logpdf(Poisson(2.0), 2) - logpdf(Poisson(4.0), 4)
    @test MeasureBase.logdensity_rel(p, AbstractMeasure(Poisson(4.0)), 3.0) == -Inf
    mn = AbstractMeasure(Multinomial(3, [0.2, 0.8]))
    X = [0.0 2.0 4.0; 6.0 4.0 2.0]
    for f in (Mul(2.0), Base.BroadcastFunction(Base.Fix1(*, 2.0)))
        @test MeasureBase.logdensities(MeasureBase.pushfwd(f, mn), X) ≈ logpdf(Multinomial(3, [0.2, 0.8]), X ./ 2)
    end
end

@testset "context precision of draws" begin
    using LinearAlgebra: I
    using StableRNGs: StableRNG
    using Statistics: mean
    using MeasureBase: GenContext, batched_rand_impl
    using Distributions: Normal, Gamma, MvNormal, Dirichlet, Poisson

    ctx(T) = GenContext{T}(StableRNG(8327))
    is_f32(x) = all(v -> Float32(v) == v, x)
    # Distributions with lower-precision parameters draw in the precision of
    # the context, not only in its type:
    for d in (Normal(0.3f0, 1.7f0), Gamma(2f0, 1.5f0), MvNormal(zeros(Float32, 2), Matrix{Float32}(I, 2, 2)), Dirichlet(Float32[1, 2, 3]))
        m = AbstractMeasure(d)
        X = batched_rand_impl(ctx(Float64), m, (200,))
        @test eltype(X) == Float64 && !is_f32(X)
        xs = [rand(ctx(Float64), m) for _ in 1:20]
        @test eltype(eltype(xs)) == Float64 && !is_f32(reduce(vcat, xs))
        @test eltype(batched_rand_impl(ctx(Float32), m, (20,))) == Float32
    end
    for d in (Normal(0.3, 1.7), Gamma(2.0, 1.5), Dirichlet([1.0, 2.0, 3.0]))
        @test eltype(batched_rand_impl(ctx(Float32), AbstractMeasure(d), (20,))) == Float32
    end
    @test mean(rand(ctx(Float64), AbstractMeasure(Gamma(2f0, 1.5f0))^4000)) ≈ 3 atol = 5 * sqrt(4.5 / 4000)
    # Without standard transports, Distributions' samplers are used:
    @test rand(ctx(Float64), AbstractMeasure(Poisson(2f0))) isa Int
end
