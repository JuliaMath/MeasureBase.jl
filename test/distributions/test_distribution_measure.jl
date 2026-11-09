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
    # Discrete measures have no transports, neither do binds of them:
    f = MeasureBase.transport_to(mbind(n -> MeasureBase.StdUniform()^n, AbstractMeasure(Poisson(3.0)), vcat), μ)
    @test_throws ArgumentError MeasureBase.transports(f, xs)
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

@testset "densities outside the support of generic distributions" begin
    using LinearAlgebra: I
    using Distributions: InverseGamma, Wishart, InverseWishart, LKJ, DirichletMultinomial, MatrixNormal, Binomial
    using MeasureBase: logdensities
    # Distributions throw outside of the support or return finite values or
    # NaN, the densities are -Inf there, point and batched:
    I2 = Matrix(1.0I, 2, 2)
    not_pd = [1.0 2.0; 2.0 1.0]
    for (d, x) in [
        (InverseGamma(2.0, 1.0), -1.0), (InverseGamma(2f0, 1f0), -1f0), (Binomial(10, 0.3), 11), (Binomial(10, 0.3), 2.5),
        (Wishart(3.0, I2), not_pd), (Wishart(3.0, I2), [1.0 0.5; 0.4 1.0]), (InverseWishart(3.0, I2), not_pd), (LKJ(2, 1.0), not_pd),
        (DirichletMultinomial(3, [1.0, 2.0]), [1.5, 1.5]), (MatrixNormal(zeros(2, 2), I2, I2), [Inf 0; 0 0]),
    ]
        m = AbstractMeasure(d)
        @test !MeasureBase.insupportof(m, x)
        @test logdensityof(m, x) == -Inf
        @test logdensityof(m, x) isa float(promote_type(eltype(x), Distributions.partype(d)))
    end
    m = AbstractMeasure(InverseGamma(2.0, 1.0))
    @test logdensities(m, [-1.0, 0.5]) ≈ [-Inf, Distributions.logpdf(InverseGamma(2.0, 1.0), 0.5)]
    @test logdensityof(m^2, [-1.0, 0.5]) == -Inf
    W = Wishart(3.0, I2)
    X = cat([2.0 0.5; 0.5 1.0], not_pd, [1.5 0.0; 0.0 1.2]; dims = 3)
    @test logdensities(AbstractMeasure(W), X) ≈ [Distributions.logpdf(W, X[:, :, 1]), -Inf, Distributions.logpdf(W, X[:, :, 3])]
    @test logdensities(AbstractMeasure(DirichletMultinomial(3, [1.0, 2.0])), [1 1.5; 2 1.5]) ≈ [Distributions.logpdf(DirichletMultinomial(3, [1.0, 2.0]), [1, 2]), -Inf]
end

# Distributions' interface doesn't require `insupport` of multivariate
# distributions:
struct _NoInsupportMvDist <: Distributions.ContinuousMultivariateDistribution end
Base.length(::_NoInsupportMvDist) = 2
Distributions._logpdf(::_NoInsupportMvDist, x::AbstractVector{<:Real}) = -sum(abs2, x) / 2 - log(2π)
Distributions._rand!(rng::MeasureBase.Random.AbstractRNG, ::_NoInsupportMvDist, x::AbstractVector{<:Real}) = MeasureBase.Random.randn!(rng, x)

@testset "densities of distributions without insupport" begin
    using MeasureBase: asmeasure, logdensities, NoFastInsupport
    d = _NoInsupportMvDist()
    m = asmeasure(d)
    x = [0.3, -0.7]
    @test MeasureBase.insupportof(m, x) isa NoFastInsupport && MeasureBase.insupportof(d, x) isa NoFastInsupport
    @test logdensityof(m, x) ≈ Distributions.logpdf(d, x)
    @test logdensityof(m, [Inf, 0.0]) == -Inf
    X = [0.3 Inf 1.2; -0.7 0.0 0.4]
    @test logdensities(m, X) ≈ [Distributions.logpdf(d, X[:, 1]), -Inf, Distributions.logpdf(d, X[:, 3])]
    @test logdensities(m, X[:, [1, 3]]) ≈ Distributions.logpdf(d, X[:, [1, 3]])
end

@testset "root measures of discrete distributions over non-integers" begin
    using Distributions: DiscreteNonParametric, Categorical, Normal
    using MeasureBase: logdensity_rel, rootmeasure, logdensities, Counting
    m = AbstractMeasure(DiscreteNonParametric([0.5, 1.5], [0.4, 0.6]))
    @test logdensity_rel(m, rootmeasure(m), 0.5) ≈ log(0.4)
    # Superpositions with continuous components and point masses see the atoms:
    s = m + AbstractMeasure(Normal())
    @test logdensityof(s, 0.5) ≈ log(0.4) && logdensityof(s, 0.7) ≈ Distributions.logpdf(Normal(), 0.7)
    @test logdensities(s, [0.5, 0.7, 1.5]) ≈ [log(0.4), Distributions.logpdf(Normal(), 0.7), log(0.6)]
    @test logdensityof(m + MeasureBase.Dirac(0.5), 1.5) ≈ log(0.6)
    # Distributions over integers count on all integers:
    @test rootmeasure(AbstractMeasure(Categorical([0.2, 0.8]))) == Counting()
    @test rootmeasure(AbstractMeasure(DiscreteNonParametric([1, 3], [0.4, 0.6]))) == Counting()
end

@testset "native draws of wrapped distributions" begin
    using Distributions: DiscreteNonParametric, Normal, Gamma, Poisson, Multinomial, DirichletMultinomial, MvNormal
    using MeasureBase: GenContext, testvalue, pushfwd
    using StableRNGs: StableRNG
    using StaticArrays: SVector
    using Static: static
    using LinearAlgebra: I
    # Atoms stay exact under lower context precisions:
    m = AbstractMeasure(DiscreteNonParametric([0.1, 0.7], [0.5, 0.5]))
    ctx32 = GenContext{Float32}(StableRNG(1))
    x = rand(ctx32, m)
    @test x isa Float64 && x in (0.1, 0.7) && logdensityof(m, x) ≈ log(0.5)
    @test logdensityof(m^5, rand(ctx32, m^5)) ≈ 5 * log(0.5)
    p2 = pushfwd(AffineMaps.Mul(2.0), m)
    @test rand(ctx32, p2) isa Float64 && logdensityof(p2, rand(ctx32, p2)) ≈ log(0.5)
    # Static powers give static arrays (without allocations where
    # Distributions' sampler doesn't allocate):
    for (d, T) in ((Normal(1.0, 2.0), Float64), (Gamma(2.0, 1.5), Float64), (Poisson(2.0), Int))
        @test @inferred(rand(StableRNG(1), AbstractMeasure(d)^static(3))) isa SVector{3,T}
    end
    alloc_of(rng, p) = (rand(rng, p); @allocated rand(rng, p))
    @test alloc_of(StableRNG(1), AbstractMeasure(Normal(1.0, 2.0))^static(3)) == 0
    @test rand(GenContext{Float32}(StableRNG(1)), AbstractMeasure(Normal(1.0, 2.0))^static(3)) isa SVector{3,Float32}
    @test rand(StableRNG(1), productmeasure((AbstractMeasure(Normal(1.0, 2.0))^static(2), StdNormal()^static(2)))) isa NTuple{2,SVector{2,Float64}}
    # Test values of multinomials (Distributions' samplers reject constant
    # draws forever for larger counts), in the support:
    for d in (Multinomial(50, [0.2, 0.3, 0.5]), Multinomial(30, [0.5, 0.5]), DirichletMultinomial(100, [1.0, 2.0]))
        x = testvalue(AbstractMeasure(d))
        @test Distributions.insupport(d, x) && sum(x) == Distributions.ntrials(d)
        @test testvalue(AbstractMeasure(d)^3) == fill(x, 3)
    end
    @test testvalue(AbstractMeasure(MvNormal([1.0, 2.0], Matrix(1.0I, 2, 2)))^2) == fill([1.0, 2.0], 2)
end

@testset "distributions of non-array variates and on the simplex" begin
    using Distributions: LKJCholesky, MvLogitNormal, Dirichlet
    using MeasureBase: rootmeasure, basemeasure, getdof, fast_dof, testvalue, logdensity_rel, NoDOF
    using LinearAlgebra: Cholesky, I
    m = AbstractMeasure(LKJCholesky(3, 1.0))
    @test rand(m) isa Cholesky && testvalue(m) isa Cholesky
    @test Distributions.insupport(LKJCholesky(3, 1.0), testvalue(m))
    @test isfinite(logdensityof(m, testvalue(m)))
    @test getdof(m) isa NoDOF && fast_dof(m) isa NoDOF
    @test_throws ArgumentError rootmeasure(m)
    @test_throws ArgumentError basemeasure(m)
    # Logit-normal densities are in the coordinates of the Dirichlet ones:
    d = MvLogitNormal([0.1, 0.2], Matrix(1.0I, 2, 2))
    ml = AbstractMeasure(d)
    x = [0.2, 0.3, 0.5]
    @test rootmeasure(ml) == rootmeasure(AbstractMeasure(Dirichlet([1.0, 2.0, 3.0])))
    @test getdof(ml) == fast_dof(ml) == 2
    @test logdensity_rel(ml, AbstractMeasure(Dirichlet([1.0, 2.0, 3.0])), x) ≈ Distributions.logpdf(d, x) - Distributions.logpdf(Dirichlet([1.0, 2.0, 3.0]), x)
end

@testset "distributions on manifolds" begin
    using Distributions: LKJ, Wishart, InverseWishart, MatrixBeta, MatrixFDist, VonMisesFisher, logpdf
    using MeasureBase: rootmeasure, getdof, fast_dof, testvalue, logdensity_rel, Lebesgue, pushfwd
    using AffineMaps: Mul
    using LinearAlgebra: I
    using StableRNGs: StableRNG
    S = [2.0 0.5; 0.5 1.0]
    # Degrees of freedom are the dimensions of the manifolds:
    for (d, n) in [
        (LKJ(3, 1.5), 3), (Wishart(5.0, S), 3), (Wishart(1.0, S), 2), (InverseWishart(5.0, S), 3),
        (MatrixBeta(3, 4.0, 5.0), 6), (MatrixFDist(4.0, 5.0, Matrix(1.0I, 3, 3)), 6), (VonMisesFisher([0.6, 0.0, 0.8], 2.0), 2),
    ]
        m = AbstractMeasure(d)
        @test getdof(m) == fast_dof(m) == n
        x = rand(StableRNG(1), d)
        @test logdensityof(m, x) ≈ logpdf(d, x)
        @test_throws ArgumentError logdensity_rel(m, Lebesgue()^size(d), x)
        d isa Union{LKJ,VonMisesFisher} && @test Distributions.insupport(d, testvalue(m))
    end
    # Distributions on the same manifold share the root:
    W, X = Wishart(5.0, S), rand(StableRNG(2), Wishart(5.0, S))
    for d2 in (Wishart(6.0, 2S), InverseWishart(5.0, S))
        @test logdensity_rel(AbstractMeasure(W), AbstractMeasure(d2), X) ≈ logpdf(W, X) - logpdf(d2, X)
    end
    @test logdensity_rel(AbstractMeasure(LKJ(3, 1.5)), AbstractMeasure(LKJ(3, 0.5)), testvalue(AbstractMeasure(LKJ(3, 1.0)))) ≈
        logpdf(LKJ(3, 1.5), testvalue(AbstractMeasure(LKJ(3, 1.0)))) - logpdf(LKJ(3, 0.5), testvalue(AbstractMeasure(LKJ(3, 1.0))))
    @test_throws ArgumentError logdensity_rel(AbstractMeasure(LKJ(2, 1.5)), AbstractMeasure(W), [1.0 0.2; 0.2 1.0])
    # Pushforwards don't apply ambient volume elements (the law of 2X is
    # Wishart(5, 4S), the ambient Jacobian would be off by log(2)):
    @test_throws ArgumentError logdensity_rel(pushfwd(Mul(2.0), AbstractMeasure(W)), AbstractMeasure(Wishart(5.0, 4S)), 2X)
end

@testset "Dirichlet draws" begin
    using Distributions: Dirichlet, Beta, ccdf
    using MeasureBase: GenContext, batched_rand_impl, testvalue
    using StableRNGs: StableRNG
    # Also for small concentrations, where all gamma variates underflow
    # together: finite draws on the simplex, with the marginal tails of
    # Beta(α₁, α₀ - α₁):
    for (T, α) in ((Float64, fill(1e-3, 3)), (Float32, fill(1f-3, 3)), (Float32, fill(1f-38, 3)), (Float64, [0.3, 2.0, 4.5]))
        m = AbstractMeasure(Dirichlet(α))
        X = batched_rand_impl(GenContext{T}(StableRNG(5)), m, (10^4,))
        @test X isa Matrix{T} && size(X) == (3, 10^4)
        @test all(isfinite, X) && all(>=(0), X) && all(s -> s ≈ 1, sum(X; dims = 1))
        a = Float64.(α)
        for t in (1e-10, 0.1, 0.5, 0.9)
            @test count(>(t), X[1, :]) / 10^4 ≈ ccdf(Beta(a[1], sum(a) - a[1]), t) atol = 0.025
        end
        @test testvalue(m) == fill(1 / 3, 3)
        @test all(isfinite, rand(StableRNG(6), T, m))
    end
    # In the precision of the context:
    x = rand(GenContext{BigFloat}(StableRNG(4)), AbstractMeasure(Dirichlet([0.5, 2.0, 3.0])))
    @test x isa Vector{BigFloat} && sum(x) ≈ 1 && any(v -> Float64(v) != v, x)
end
