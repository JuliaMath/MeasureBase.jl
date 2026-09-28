using Test

using MeasureBase
using MeasureBase: pushfwd, StdUniform, StdExponential, StdLogistic
using MeasureBase: pushfwd, PushforwardMeasure
using MeasureBase: transport_to, unsafe_logdensityof
using MeasureBase: productmeasure, mbind
import Zygote
using InverseFunctions: inverse
using ChangesOfVariables: with_logabsdet_jacobian
import Statistics: var
using DensityInterface: logdensityof
using LogExpFunctions
using SpecialFunctions: erfc, erfcinv
import InverseFunctions: inverse, FunctionWithInverse, setinverse
using IrrationalConstants: invsqrt2, sqrt2
import ChangesOfVariables: with_logabsdet_jacobian
using MeasureBase.Interface: transport_to, test_transport

Φ(z) = erfc(-z * invsqrt2) / 2
Φinv(p) = -erfcinv(2 * p) * sqrt2

with_logabsdet_jacobian(f::FunctionWithInverse, x) = with_logabsdet_jacobian(f.f, x)

with_logabsdet_jacobian(::typeof(Φ), z) = (Φ(z), logdensityof(StdNormal(), z))

function with_logabsdet_jacobian(::typeof(Φinv), p)
    z = Φinv(p)
    (z, -logdensityof(StdNormal(), z))
end

inverse(::typeof(Φ)) = Φinv
inverse(::typeof(Φinv)) = Φ

var(::StdNormal) = 1.0
var(::StdExponential) = 1.0
var(::StdUniform) = 1 / 12
var(::StdLogistic) = π^2 / 3

function test_pushfwd(f, μ, ν_ref)
    @testset "pushfwd($f, $μ)" begin
        @inferred(pushfwd(f, μ))
        ν = pushfwd(f, μ)

        test_transport(ν, ν_ref)
        test_transport(ν_ref, ν)

        y = rand(ν_ref)
        @test isapprox(@inferred(logdensityof(ν, y)), logdensityof(ν_ref, y), atol = 1e-10)
        @test isapprox(
            @inferred(unsafe_logdensityof(ν, y)),
            unsafe_logdensityof(ν_ref, y),
            atol = 1e-10,
        )

        @test isapprox(var(rand(ν^(10^5))), var(ν_ref), rtol = 0.05)

        @test abs(transport_to(StdLogistic(), ν)(y)) ≈
              abs(transport_to(StdLogistic(), ν_ref)(y))
    end
end

@testset "transformedmeasure.jl" begin
    # (f, μ, ν_ref), so that pushfwd(f, μ) ≅ ν_ref
    triples = [
        ((-) ∘ log, StdUniform(), StdExponential())
        (exp ∘ (-), StdExponential(), StdUniform())
        (logit, StdUniform(), StdLogistic())
        (logistic, StdLogistic(), StdUniform())
        (Φ, StdNormal(), StdUniform())
        (Φinv, StdUniform(), StdNormal())
    ]

    for (f, μ, ν_ref) in triples
        test_pushfwd(f, μ, ν_ref)
        test_pushfwd(setinverse(f, inverse(f)), μ, ν_ref)
    end

    @testset "Pushforward-of-pushforward" begin
        for (f, μ, ν_ref) in triples
            finv = inverse(f)
            test_pushfwd(finv, pushfwd(f, μ), μ)
        end
    end
end

using IrrationalConstants: loghalf

@testset "Half" begin
    μ = Half(StdNormal())

    # Test logdensityof for positive values
    @test logdensityof(μ, 1.0) ≈ logdensityof(StdNormal(), 1.0) - loghalf

    # Test logdensityof for negative values (should return -Inf)
    @test logdensityof(μ, -1.0) == -Inf

    # Test logdensityof for zero
    @test logdensityof(μ, 0.0) ≈ logdensityof(StdNormal(), 0.0) - loghalf

    @test logdensities(μ, [-1.0, 1.0]) ≈ [-Inf, logdensityof(μ, 1.0)]
    @test MeasureBase.insupport(μ, -1.0) == false && MeasureBase.insupport(μ, 1.0) == true
    @test logdensityof(Half(2 * StdNormal()), 1.0) ≈ log(2) + logdensityof(μ, 1.0)

    # Half requires scalar variates and parents symmetric about zero:
    @test_throws ArgumentError Half(StdExponential())
    @test_throws ArgumentError Half(StdUniform())
    @test_throws ArgumentError Half(StdNormal()^3)
end

using Test
using MeasureBase
using MeasureBase: gettransform
using InverseFunctions
using StaticArrays
using ChangesOfVariables
using ChangesOfVariables: setladj

@testset "PushforwardMeasure" begin
    # Test basic pushforward construction
    μ = StdNormal()
    f = exp
    ν = pushfwd(f, μ)

    @test ν isa PushforwardMeasure

    # TODO: How do we access the parent?
    # @test parent(ν) === μ
    @test gettransform(ν) === f

    # Test logdensity with AdaptRootMeasure
    x = 0.5
    y = f(x)
    ld = logdensityof(ν, y)
    @test !isnan(ld)
    @test isfinite(ld)

    # Test logdensity with PushfwdRootMeasure
    ν_no_corr = pushfwd(f, μ, PushfwdRootMeasure())
    ld_no_corr = logdensityof(ν_no_corr, y)
    @test !isnan(ld_no_corr)
    @test isfinite(ld_no_corr)

    # Test non-bijective pushforward
    g = x -> x^2
    @test MeasureBase.getdof(pushfwd(g, μ)) isa MeasureBase.NoDOF
    @test_throws ArgumentError logdensityof(pushfwd(g, μ), 1.0)
    @test_throws ArgumentError logdensityof(pushfwd(g, μ, PushfwdRootMeasure()), 1.0)

    # Test edge cases for _combine_logd_with_ladj
    @test MeasureBase._combine_logd_with_ladj(-Inf, Inf) == -Inf  # Zero density wins
    @test MeasureBase._combine_logd_with_ladj(1.0, -Inf) == -Inf
    @test logdensityof(pushfwd(log, StdExponential()), -Inf) == -Inf

    # Test composition of pushforwards
    h = x -> exp(x)
    ν_comp = pushfwd(h, ν)
    @test ν_comp isa PushforwardMeasure

    # TODO: How do we access the parent?
    # @test parent(ν_comp) === μ

    # Test identity pushforward
    ν_id = pushfwd(identity, μ)
    @test ν_id === μ

    # Test rootmeasure behavior
    @test rootmeasure(ν) === rootmeasure(μ)  # AdaptRootMeasure
    @test rootmeasure(ν_no_corr) isa PushforwardMeasure  # PushfwdRootMeasure

    # Pushforwards adapting the root have the origin's root as their base,
    # pushforwards of the root pull weights out of the pushed-forward base:
    @test basemeasure(ν) === rootmeasure(μ)
    @test basemeasure(ν_no_corr) isa WeightedMeasure{<:Any,<:PushforwardMeasure}
    @test basemeasure(ν_no_corr).base.style isa PushfwdRootMeasure

    # Both styles describe the same measure:
    @test logdensity_rel(ν, ν_no_corr, 2.0) ≈ 0 atol = 1e-12
    @test logdensity_rel(ν_no_corr, ν, 2.0) ≈ 0 atol = 1e-12
    @test logdensity_rel(ν, rootmeasure(μ), 2.0) ≈ logdensityof(μ, log(2.0)) - log(2.0)
    E = Base.BroadcastFunction(exp)
    νv, νv_no_corr = pushfwd(E, μ^2), pushfwd(E, μ^2, PushfwdRootMeasure())
    @test logdensity_rel(νv, νv_no_corr, [2.0, 3.0]) ≈ 0 atol = 1e-12
    @test logdensities(νv, [2.0 1.0; 3.0 0.5]) ≈ [logdensityof(νv, [2.0, 3.0]), logdensityof(νv, [1.0, 0.5])]

    # Rank-changing pushforwards reshape the root:
    mr = MeasureBase.mreshape(μ^6, (2, 3))
    @test basemeasure(mr) == MeasureBase.LebesgueBase()^(2, 3)
    xr = reshape(collect(0.1:0.1:0.6), 2, 3)
    @test logdensityof(mr, xr) ≈ logdensityof(μ^6, vec(xr))

    # Mixtures of pushforward location families:
    shifted = pushfwd(setinverse(setladj(x -> x + 1, zero), setladj(y -> y - 1, zero)), μ)
    s = superpose([μ, shifted])
    @test logdensityof(s, 0.3) ≈ log(exp(logdensityof(μ, 0.3)) + exp(logdensityof(μ, -0.7)))
    @test sum(x -> exp(logdensityof(s, x)), -12:0.01:13) * 0.01 ≈ 2

    # Test massof
    # TODO: mass interface is very incomplete
    # @test massof(ν) == massof(μ)

    # Test rand
    @test rand(ν) isa Real
    @test insupport(ν, rand(ν)) != false

    # Test pullback
    pb = pullbck(f, ν)
    @test pb isa PushforwardMeasure
    @test logdensityof(pb, y) ≈ logdensityof(μ, y)

    # Test deprecated pullback
    @test_deprecated pullback(f, μ)
end

@testset "PushFwdStyle types" begin
    @test AdaptRootMeasure() isa PushFwdStyle
    @test PushfwdRootMeasure() isa PushFwdStyle
    @test MeasureBase.WithVolCorr === AdaptRootMeasure
    @test MeasureBase.NoVolCorr === PushfwdRootMeasure

    @testset "output size of pushforwards of tuple products" begin
        Pt = productmeasure((StdNormal(), StdUniform()^2))
        ν = pushfwd(x -> vcat(x[1], x[2]), Pt)
        @test MeasureBase.mspace_flatsize(ν) == (3,)
        @test MeasureBase.mspace_ndims(typeof(ν)) == 1
        @test MeasureBase.mspace_flatsize(pushfwd(x -> x, mbind(x -> StdNormal()^(x > 0 ? 1 : 2), StdNormal()))) isa MeasureBase.NoMSpaceElementSize
    end

    @testset "construction inside differentiated functions" begin
        flat(x::Tuple) = vcat(x[1], x[2])
        unflat(v::AbstractVector) = (v[1], v[2:end])
        InverseFunctions.inverse(::typeof(flat)) = unflat
        InverseFunctions.inverse(::typeof(unflat)) = flat
        ChangesOfVariables.with_logabsdet_jacobian(::typeof(flat), x) = (flat(x), zero(eltype(x[2])))
        ChangesOfVariables.with_logabsdet_jacobian(::typeof(unflat), v) = (unflat(v), zero(eltype(v)))
        Pt = productmeasure((StdNormal(), StdUniform()^2))
        g(v) = logdensityof(pushfwd(flat, Pt), v)
        v = [0.3, 0.2, 0.7]
        @test Zygote.gradient(g, v)[1] ≈ [-0.3, 0.0, 0.0]
    end
end
