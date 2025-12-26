using Test

using MeasureBase
using MeasureBase: pushfwd, StdUniform, StdExponential, StdLogistic
using MeasureBase: pushfwd, PushforwardMeasure
using MeasureBase: transport_to, unsafe_logdensityof
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
        x = inverse(f)(y)

        β = basemeasure(μ)
        @test isapprox(
            logdensity_rel(ν, pushfwd(f, β), y),
            logdensity_rel(μ, β, x),
            atol = 1e-10,
        )

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
end

using Test
using MeasureBase
using MeasureBase: gettransform
using InverseFunctions
using StaticArrays
using ChangesOfVariables

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
    @test_throws ArgumentError logdensityof(pushfwd(g, μ), 1.0)
    @test_throws ArgumentError logdensityof(pushfwd(g, μ, PushfwdRootMeasure()), 1.0)

    # Test edge cases for _combine_logd_with_ladj
    @test MeasureBase._combine_logd_with_ladj(-Inf, Inf) == -Inf  # Zero density wins
    @test MeasureBase._combine_logd_with_ladj(1.0, -Inf) ≈ MeasureBase.near_neg_inf(Float64)

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

    # Test basemeasure
    @test basemeasure(ν) isa PushforwardMeasure
    @test basemeasure(ν).style isa PushfwdRootMeasure

    # Test massof
    # TODO: mass interface is very incomplete
    # @test massof(ν) == massof(μ)

    # Test rand
    @test rand(ν) isa Real
    @test insupport(ν, rand(ν))

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
end
