using Test

using MeasureBase: pushfwd, StdUniform, StdExponential, StdLogistic
using MeasureBase: pushfwd, PushforwardMeasure
using MeasureBase: transport_to, unsafe_logdensityof
import Statistics: var
using DensityInterface: logdensityof
using LogExpFunctions
using SpecialFunctions: erfc, erfcinv
import InverseFunctions: inverse


Φ(z) = erfc(-z * invsqrt2)/2
Φinv(p) = -erfcinv(2*p) * sqrt2

inverse(::typeof(Φ)) = Φinv
inverse(::typeof(Φinv)) = Φ

var(::StdNormal) = 1.0
var(::StdExponential) = 1.0
var(::StdUniform) = 1/12
var(::StdLogistic) = π^2/3

@testset "transformedmeasure.jl" begin
    for (f, μ, ν_ref) in [
        ((-) ∘ log1p ∘ (-), StdUniform(), StdExponential())
        ((-) ∘ log, StdUniform(), StdExponential())
        (exp ∘ (-), StdExponential(), StdUniform())
        (logit, StdUniform(), StdLogistic())
        (logistic, StdLogistic(), StdUniform())
        # (Φ, StdNormal(), StdUniform())
        # (Φinv, StdUniform(), StdNormal())
        ]

        finv = inverse(f)

        @testset "pushfwd($f, $μ)" begin
            @test @inferred(pushfwd(f, μ)) isa PushforwardMeasure
            ν = pushfwd(f, μ)

            y = rand(ν_ref)
            @test @inferred(logdensityof(ν, y)) ≈ logdensityof(ν_ref, y)
            @test @inferred(unsafe_logdensityof(ν, y)) ≈ unsafe_logdensityof(ν_ref, y)

            @test isapprox(var(rand(ν^(10^5))), var(ν_ref), rtol = 0.05)

            @test abs(transport_to(StdLogistic(), ν)(y)) ≈ abs(transport_to(StdLogistic(), ν_ref)(y))
        end

        @testset "pushfwd($finv, pushfwd($f, μ))" begin
            # TODO
        end
    end
end
