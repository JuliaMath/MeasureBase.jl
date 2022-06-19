using Test

using MeasureBase: pushfwd, StdUniform, StdExponential, StdLogistic
using MeasureBase: pushfwd, PushforwardMeasure
using MeasureBase: transport_to
using Statistics: var
using DensityInterface: logdensityof

@testset "transformedmeasure.jl" begin
    μ = StdUniform()
    @test @inferred(pushfwd((-) ∘ log1p ∘ (-), μ)) isa PushforwardMeasure
    ν = pushfwd((-) ∘ log1p ∘ (-), μ)
    ν_ref = StdExponential()

    y = rand(ν_ref)
    @test @inferred(logdensityof(ν, y)) ≈ logdensityof(ν_ref, y)
    
    @test isapprox(var(rand(ν^(10^5))), 1, rtol = 0.05)

    @test transport_to(StdLogistic(), ν)(y) ≈ transport_to(StdLogistic(), ν)(y)
end
