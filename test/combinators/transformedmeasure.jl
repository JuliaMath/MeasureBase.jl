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
