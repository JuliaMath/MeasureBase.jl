using Test
using MeasureBase

@testset "public API" begin
    @static if VERSION >= v"1.11"
        public_names = filter(n -> Base.ispublic(MeasureBase, n), names(MeasureBase; all = true))
        @test all(n -> isdefined(MeasureBase, n), public_names)
        # Extension points and internals used by downstream packages:
        for n in (
            :logdensityof_impl, :batched_logdensityof_impl, :rand_impl, :batched_rand_impl,
            :transport_to_std, :transport_from_std, :preferred_stdmeasure, :mspace_ndims,
            :mspace_flatsize, :checked_arg, :NoDOF, :UnknownMass, :DensityMeasure,
            :TransportFunction, :GenContext,
        )
            @test Base.ispublic(MeasureBase, n) && !Base.isexported(MeasureBase, n)
        end
    end
end

@testset "conditioning" begin
    P = productmeasure((a = StdNormal(), b = StdExponential()))
    C = P | (b = 1.0,)
    @test logdensityof(C, (a = 0.3, b = 1.0)) ≈ logdensityof(StdNormal(), 0.3)
    @test logdensityof(C, (a = 0.3, b = 2.0)) == -Inf
    @test rand(C).b == 1.0
    @test (StdNormal() | isfinite) == MeasureBase.condition(StdNormal(), isfinite)
end

@testset "likelihoods" begin
    k(σ) = weightedmeasure(-log(σ), StdNormal())
    ℓ = likelihoodof(k, 0.5)
    @test logdensityof(ℓ, 2.0) ≈ logdensityof(k(2.0), 0.5)
    @test likelihood_kernel(ℓ) === k && likelihood_obs(ℓ) == 0.5
    for l in (
        Base.Fix2(densityof, 0.5) ∘ k,
        MeasureBase.funcdensity(Base.Fix2(densityof, 0.5) ∘ k),
        MeasureBase.logfuncdensity(Base.Fix2(logdensityof, 0.5) ∘ k),
    )
        @test as_likelihood(l) == ℓ
    end
    @test as_likelihood(ℓ) === ℓ

    # A likelihood as integrand gives an unnormalized posterior:
    post = mintegrate(ℓ, StdExponential())
    @test logdensityof(post, 2.0) ≈ logdensityof(ℓ, 2.0) + logdensityof(StdExponential(), 2.0)
end

@testset "mintegrate" begin
    m = mintegrate(x -> 2x, StdExponential())
    @test logdensityof(m, 1.5) ≈ log(3.0) - 1.5
    @test logdensityof(m, -1.5) == -Inf
    @test logdensityof(∫(x -> x^2, Lebesgue()), 3.0) ≈ log(9.0)
    @test logdensityof(∫exp(x -> -x^2, Lebesgue()), 3.0) ≈ -9.0
end

@testset "half" begin
    h = half(StdNormal())
    @test logdensityof(h, 0.5) ≈ log(2) + logdensityof(StdNormal(), 0.5)
    @test logdensityof(h, -0.5) == -Inf
    @test all(≥(0), [rand(h) for _ in 1:20])
end

@testset "dirac" begin
    @test dirac(Dirac(2.0)) == Dirac(2.0)
    d = dirac(StdNormal())
    @test d isa Dirac && logdensityof(d, d.x) == 0
end

@testset "valdomain" begin
    @test valdomain(1.0) isa MeasureBase.UnknownDomain
end

@testset "unsafe_logdensity_rel" begin
    μ, ν = StdNormal(), StdUniform()
    @test unsafe_logdensity_rel(μ, ν, 0.3) ≈ logdensity_rel(μ, ν, 0.3)
    # No support checks, unlike logdensity_rel:
    @test logdensity_rel(μ, ν, 1.5) == Inf
    @test unsafe_logdensity_rel(μ, ν, 1.5) ≈ logdensityof(μ, 1.5)
end

@testset "hash follows equality" begin
    makers = (
        () -> StdNormal()^(2, 3),
        () -> weightedmeasure(0.5, StdNormal()),
        () -> Dirac([1.0, 2.0]),
        () -> superpose(StdNormal(), Dirac(0.0)),
        () -> pushfwd(exp, StdNormal()),
        () -> mintegrate_exp(sin, StdNormal()),
        () -> restrict(>(0), StdNormal()),
        () -> half(StdNormal()),
        () -> mbind(Dirac, StdNormal()),
        () -> mcombine(vcat, StdNormal()^2, StdExponential()^3),
    )
    for mk in makers
        a, b = mk(), mk()
        @test a == b && hash(a) == hash(b)
    end
    @test Dirac([1.0, 2.0]) != Dirac([1.0, 3.0])
    @test weightedmeasure(0.5, StdNormal()) != weightedmeasure(0.6, StdNormal())
end
