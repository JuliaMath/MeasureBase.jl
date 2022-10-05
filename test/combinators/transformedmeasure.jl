using Test

using MeasureBase: pushfwd, StdUniform, StdExponential, StdLogistic
using MeasureBase: pushfwd, PushforwardMeasure
using MeasureBase: transport_to, unsafe_logdensityof
import Statistics: var
using DensityInterface: logdensityof
using LogExpFunctions

var(::StdNormal) = 1.0
var(::StdExponential) = 1.0
var(::StdUniform) = 1/12

@testset "transformedmeasure.jl" begin
    for (f, μ, ν_ref) in [
        ((-) ∘ log1p ∘ (-), StdUniform(), StdExponential())
        ((-) ∘ log, StdUniform(), StdExponential())
        (logit, StdUniform(), StdLogistic())
        (logistic, StdLogistic(), StdUniform())
        ]

        @testset "pushfwd($f, $μ)" begin
            @test @inferred(pushfwd(f, μ)) isa PushforwardMeasure
            ν = pushfwd(f, μ)

            y = rand(ν_ref)
            @test @inferred(logdensityof(ν, y)) ≈ logdensityof(ν_ref, y)
            @test @inferred(unsafe_logdensityof(ν, y)) ≈ unsafe_logdensityof(ν_ref, y)

            @test isapprox(var(rand(ν^(10^5))), var(ν_ref), rtol = 0.05)

            @test transport_to(StdLogistic(), ν)(y) ≈ transport_to(StdLogistic(), ν)(y)
        end
    end
end
