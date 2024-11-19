# test/combinators/conditional.jl
using Test
using MeasureBase
using Random: MersenneTwister

@testset "Conditional" begin
    # Create a simple conditional measure
    base_measure = StdNormal()
    condition(x) = abs(x) <= 2  # Only accept values in [-2, 2]
    
    cond_measure = @inferred Conditional(base_measure, condition)
    
    # Test basic properties
    @test basemeasure(cond_measure) === base_measure
    
    # Test sampling with rejection sampling
    rng = MersenneTwister(123)
    samples = [rand(rng, cond_measure) for _ in 1:100]
    @test all(condition, samples)
    
    # Test density
    x = 1.0
    @test logdensityof(cond_measure, x) ≈ logdensityof(base_measure, x)
    @test logdensityof(cond_measure, 3.0) == -Inf  # Outside condition
    
    # Test support
    @test insupport(cond_measure, 1.0)
    @test !insupport(cond_measure, 3.0)
end