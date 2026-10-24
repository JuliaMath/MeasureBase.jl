# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

import Distributions
using Distributions: Distribution
import MeasureBase
import AffineMaps
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
end

@testset "Dirichlet reference measure" begin
    using Distributions: Dirichlet, MvNormal, logpdf
    using LinearAlgebra: I
    d1, d2 = AbstractMeasure(Dirichlet(ones(3))), AbstractMeasure(Dirichlet([2.0, 3.0, 4.0]))
    L = MeasureBase.Lebesgue(MeasureBase.StdSimplex(3))
    @test MeasureBase.basemeasure(d1) == L && MeasureBase.rootmeasure(d1) == L
    @test MeasureBase.getdof(d1) == MeasureBase.getdof(L) == 2
    @test MeasureBase.fast_dof(d1) == 2
    x = [0.2, 0.3, 0.5]
    @test MeasureBase.logdensity_rel(d2, d1, x) ≈ logpdf(Dirichlet([2.0, 3.0, 4.0]), x) - log(2)
    @test MeasureBase.logdensity_rel(d2, L, x) ≈ logpdf(Dirichlet([2.0, 3.0, 4.0]), x)
    # Singular relative to measures with densities on ℝ³:
    for ν in (AbstractMeasure(MvNormal(zeros(3), I)), MeasureBase.Lebesgue()^3)
        @test_throws ArgumentError MeasureBase.logdensity_rel(d1, ν, x)
        @test_throws ArgumentError MeasureBase.logdensity_rel(ν, d1, x)
    end

    # Maps of the ambient space push the simplex reference forward, they
    # have no volume element for it:
    p = MeasureBase.pushfwd(AffineMaps.Mul(2I), d2)
    @test MeasureBase.rootmeasure(p) == MeasureBase.pushfwd(AffineMaps.Mul(2I), L)
    @test MeasureBase.logdensityof(p, 2x) ≈ logpdf(Dirichlet([2.0, 3.0, 4.0]), x)
    @test MeasureBase.logdensityof(p, 3x) == -Inf
    @test MeasureBase.logdensity_rel(p, MeasureBase.pushfwd(AffineMaps.Mul(2I), d1), 2x) ≈ logpdf(Dirichlet([2.0, 3.0, 4.0]), x) - log(2)
    @test sum(rand(p)) ≈ 2
    # Pushforwards by transports have the root measure of the target:
    U2 = MeasureBase.StdUniform()^2
    q = MeasureBase.pushfwd(MeasureBase.transport_to(d2, U2), U2)
    @test MeasureBase.rootmeasure(q) == L
    @test MeasureBase.logdensityof(q, x) ≈ logpdf(Dirichlet([2.0, 3.0, 4.0]), x)
    @test MeasureBase.logdensity_rel(q, d2, x) ≈ 0 atol = 1e-12
    @test MeasureBase.logdensity_rel(q, d1, x) ≈ logpdf(Dirichlet([2.0, 3.0, 4.0]), x) - log(2)
    r = MeasureBase.pushfwd(MeasureBase.transport_to(U2, d2), d2)
    u = MeasureBase.transport_to(U2, d2)(x)
    @test MeasureBase.rootmeasure(r) == MeasureBase.rootmeasure(U2)
    @test MeasureBase.logdensity_rel(r, U2, u) ≈ 0 atol = 1e-12
end

@testset "fast DOF of wrapped distributions" begin
    using Distributions: Normal, Poisson, Dirichlet, MvNormal
    using MeasureBase: fast_dof, restrict
    @test @inferred(fast_dof(AbstractMeasure(Normal(0.3, 2.0)))) === MeasureBase.static(1)
    @test @inferred(fast_dof(AbstractMeasure(Poisson(2.0)))) === MeasureBase.static(1)
    @test fast_dof(restrict(x -> x > 1, AbstractMeasure(Poisson(2.0)))) === MeasureBase.static(1)
    @test fast_dof(AbstractMeasure(Dirichlet([1.0, 2.0, 3.0, 4.0]))) == 3
    @test fast_dof(AbstractMeasure(MvNormal([1.0 0.2; 0.2 2.0]))) == 2
    @test fast_dof(AbstractMeasure(Normal())^3) == 3
    # Transport partners take the fast DOF, static where it is static:
    @test transport_to(MeasureBase.StdNormal, Normal(0.3, 2.0)).ν === MeasureBase.StdNormal()
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
    # Origins without random variates:
    pr = MeasureBase.pushfwd(Mul(2.0), MeasureBase.restrict(x -> x > 1, m))
    @test MeasureBase.logdensityof(pr, 4.0) == MeasureBase.logdensityof(m, 2.0)
    @test MeasureBase.logdensityof(pr, 2.0) == -Inf
    @test MeasureBase.logdensities(pr, [2.0, 4.0, 6.0]) ≈ [-Inf, logpdf(Poisson(2.0), 2), logpdf(Poisson(2.0), 3)]
    mn = AbstractMeasure(Multinomial(3, [0.2, 0.8]))
    X = [0.0 2.0 4.0; 6.0 4.0 2.0]
    for f in (Mul(2.0), Base.BroadcastFunction(Base.Fix1(*, 2.0)))
        @test MeasureBase.logdensities(MeasureBase.pushfwd(f, mn), X) ≈ logpdf(Multinomial(3, [0.2, 0.8]), X ./ 2)
    end
end

@testset "Poisson process as a bind" begin
    using Distributions: Poisson, Normal, logpdf
    using MeasureBase: StdNormal, mbind
    using Random: Xoshiro
    # Poisson counts are floats, powers take integer-valued counts:
    μ = mbind(n -> StdNormal()^n, AbstractMeasure(Poisson(3.0)), vcat)
    ℓ_ref(x) = logpdf(Poisson(3.0), x[1]) + sum(logpdf.(Normal(), x[2:end]); init = 0.0)
    xs = [rand(Xoshiro(i), μ) for i in 1:20]
    @test all(x -> length(x) == x[1] + 1, xs) && any(x -> x[1] >= 3, xs)
    @test rand(Xoshiro(1), mbind(n -> StdNormal()^n, AbstractMeasure(Poisson(1e-3)), vcat)) == [0.0]
    @test MeasureBase.logdensityof(μ, xs[1]) ≈ ℓ_ref(xs[1])
    @test MeasureBase.logdensityof(μ, [0.0]) ≈ -3
    push!(xs, [0.0])
    @test MeasureBase.logdensities(μ, xs) ≈ ℓ_ref.(xs)
    @test MeasureBase.logdensities(μ, [2.0 2.0; 0.1 -0.3; 0.2 0.4]) ≈ ℓ_ref.([[2.0, 0.1, 0.2], [2.0, -0.3, 0.4]])
    @test MeasureBase.logdensityof(μ, [2.5, 0.1, 0.2]) == -Inf
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

@testset "test values" begin
    # Constant RNGs sample univariate distributions by inversion, at the
    # median, in the type of their native draws:
    for d in (DiscreteUniform(1, 6), Poisson(3.5), Bernoulli(0.3), Categorical([0.2, 0.5, 0.3]), Gamma(2.0, 3.0))
        m = AbstractMeasure(d)
        x = MeasureBase.testvalue(m)
        @test x isa eltype(d) && x == quantile(d, 0.5)
        @test Distributions.insupport(d, x) && isfinite(logdensityof(m, x))
        @test MeasureBase.testvalue(Float32, m) ≈ x
        @test rand(MeasureBase.ConstantRNG(), Float64, m^3) == fill(x, 3)
    end
end
