using Test

using MeasureBase: AbstractMeasure
using MeasureBase: StdExponential, StdLogistic, StdNormal, StdUniform
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

module _TestOperatorExports
    using Test
    using MeasureBase

    @testset "operators via using MeasureBase" begin
        μ, ν = StdNormal(), StdExponential()
        @test ⊗(μ, ν) == productmeasure((μ, ν))
        @test logdensityof(∫exp(x -> -x, μ), 0.5) ≈ logdensityof(μ, 0.5) - 0.5
        @test logdensityof(∫(x -> 2.0, μ), 0.5) ≈ logdensityof(μ, 0.5) + log(2)
        @test log𝒹(μ, Lebesgue())(0.5) ≈ logdensityof(μ, 0.5)
        @test 𝒹(μ, Lebesgue())(0.5) ≈ exp(logdensityof(μ, 0.5))
        @test (exp ⋄ μ) == pushfwd(exp, μ)
        @test (μ ⊙ exp) == pullbck(exp, μ)
        @test (μ ▷ Dirac) == mbind(Dirac, μ)
    end
end
