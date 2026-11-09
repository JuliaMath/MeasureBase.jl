using Test
using Random
using MeasureBase
using MeasureBase.Interface: test_interface, test_transport
import Distributions

@testset "public API" begin
    @static if VERSION >= v"1.11"
        public_names = filter(n -> Base.ispublic(MeasureBase, n), names(MeasureBase; all = true))
        @test all(n -> isdefined(MeasureBase, n), public_names)
        # Extension points and internals used by downstream packages:
        for n in (
            :logdensityof_impl, :batched_logdensityof_impl, :rand_impl, :batched_rand_impl,
            :transport_to_std, :transport_from_std, :preferred_stdmeasure, :mspace_ndims,
            :mspace_flatsize, :stream_ndims, :SizeLike, :size_tuple, :static, :checked_arg, :NoDOF,
            :UnknownMass, :DensityMeasure, :TransportFunction, :GenContext, :Bind, :MKernel,
            :kernelmeasure, :valuecombinator, :combinedmeasures, :Interface, :convert_realtype,
        )
            @test Base.ispublic(MeasureBase, n) && !Base.isexported(MeasureBase, n)
        end
    end
end

@testset "deprecated names" begin
    # Deprecated names are reachable, but not exported:
    for n in (:insupport, :params, :pullback, :density_def)
        @test isdefined(MeasureBase, n) && !Base.isexported(MeasureBase, n)
    end
    @test (@test_deprecated MeasureBase.insupport(StdUniform(), 1.5)) == false
    @test (@test_deprecated MeasureBase.insupport(StdUniform())(0.5)) == true
    @test (@test_deprecated MeasureBase.params(StdNormal())) == NamedTuple()
    @test (@test_deprecated MeasureBase.pullback(exp, StdNormal())) == pullbck(exp, StdNormal())
    @test (@test_deprecated MeasureBase.density_def(StdNormal(), 0.3)) ≈ exp(logdensity_def(StdNormal(), 0.3))
    # so that they don't clash with the names of Distributions:
    m = Module()
    @eval m using MeasureBase, Distributions
    @test (@eval m insupport) === Distributions.insupport
    @test (@eval m params) === Distributions.params
end

@testset "gentype" begin
    # MeasureBase doesn't extend `Random.gentype`, the type of `rand` depends
    # on the generative context:
    @test !Base.isexported(MeasureBase, :gentype)
    @test !any(m -> parentmodule(m) === MeasureBase, methods(Random.gentype))
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


# Records test results without reporting them, to test the test helpers:
struct _RecordingTestSet <: Test.AbstractTestSet
    description::String
    results::Vector{Any}
end
_RecordingTestSet(description::AbstractString; kwargs...) = _RecordingTestSet(description, [])
Test.record(ts::_RecordingTestSet, result) = (push!(ts.results, result); result)
Test.record(ts::_RecordingTestSet, child::_RecordingTestSet) = (append!(ts.results, child.results); child)
function Test.finish(ts::_RecordingTestSet)
    parent = Test.get_testset_depth() > 0 ? Test.get_testset() : nothing
    parent isa _RecordingTestSet && Test.record(parent, ts)
    return ts
end

# A normal distribution whose batched density and batched transport to
# StdNormal disagree with the point kernels:
struct _InconsistentNormal <: AbstractMeasure end
MeasureBase.logdensity_def(::_InconsistentNormal, x) = -x^2 / 2
MeasureBase.basemeasure(::_InconsistentNormal) = weightedmeasure(-log(2π) / 2, Lebesgue())
MeasureBase.insupportof(::_InconsistentNormal, x) = true
MeasureBase.massof(::_InconsistentNormal) = MeasureBase.static(1.0)
MeasureBase.mspace_flatsize(::Type{<:_InconsistentNormal}) = ()
MeasureBase.mspace_flatsize(::_InconsistentNormal) = ()
MeasureBase.preferred_stdmeasure(::Type{<:_InconsistentNormal}) = StdNormal
MeasureBase.transport_to_std(::Type{StdNormal}, ::_InconsistentNormal, x) = x
MeasureBase.transport_from_std(::Type{StdNormal}, ::_InconsistentNormal, z) = z
MeasureBase.batched_logdensityof_impl(::_InconsistentNormal, X) = @. -X^2 - log(2π) / 2
MeasureBase.batched_transport_to_std(::Type{StdNormal}, ::_InconsistentNormal, X) = 2 .* reshape(X, 1, size(X)...)

@testset "Interface test helpers check batched kernels" begin
    for check in (() -> test_interface(_InconsistentNormal()), () -> test_transport(StdNormal(), _InconsistentNormal()))
        ts = @testset _RecordingTestSet "inconsistent batched kernels" begin
            check()
        end
        @test any(r -> r isa Test.Fail, ts.results)
        @test !any(r -> r isa Test.Error, ts.results)
    end
end
