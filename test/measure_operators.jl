using Test

using MeasureBase: AbstractMeasure
using MeasureBase: StdExponential, StdLogistic, StdUniform
using MeasureBase: pushfwd, pullbck, mbind, productmeasure
using MeasureBase: mintegrate, mintegrate_exp, density_rel, logdensity_rel
using MeasureBase.MeasureOperators: ⋄, ⊙, ▷, ⊗, ∫, ∫exp, 𝒹, log𝒹

@testset "MeasureOperators" begin
    μ = StdExponential()
    ν = StdUniform()
    k(σ) = pushfwd(x -> σ * x, StdNormal())
    μs = (StdExponential(), StdLogistic(), StdUniform())
    f = sqrt

    @test @inferred(f ⋄ μ) == pushfwd(f, μ)
    @test @inferred(ν ⊙ f) == pullbck(f, ν)
    @test @inferred(μ ▷ k) == mbind(k, μ)
    @test @inferred(⊗(μs...)) == productmeasure(μs)
    @test @inferred(∫(f, μ)) == mintegrate(f, μ)
    @test @inferred(∫exp(f, μ)) == mintegrate_exp(f, μ)
    @test @inferred(𝒹(ν, μ)) == density_rel(ν, μ)
    @test @inferred(log𝒹(ν, μ)) == logdensity_rel(ν, μ)
end
