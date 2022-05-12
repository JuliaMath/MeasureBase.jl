using Random: MersenneTwister
using Test

using MeasureBase
using MeasureBase: _logweight, weightedmeasure, WeightedMeasure

@testset "weighted" begin
    @test iszero(_logweight(Lebesgue(ℝ)))
    μ₀ = Dirac(0.0)
    w = 2.0
    μ = @inferred w * μ₀
    @test μ == WeightedMeasure(log(w), μ₀) == weightedmeasure(log(w), μ₀)
    @test μ isa WeightedMeasure
    @test _logweight(μ) == log(w)
    @test _logweight(w * μ) == 2 * log(w)
    @test rand(MersenneTwister(123), μ) == rand(MersenneTwister(123), μ₀)
    x = rand()
    @test logdensity_def(μ, x) == log(w)
    @test logdensityof(μ, x) == logdensityof(μ₀, x)
end
