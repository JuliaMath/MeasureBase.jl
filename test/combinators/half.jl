using Test
using MeasureBase
using Random: MersenneTwister
using LogExpFunctions: loghalf
using IrrationalConstants: log2π

@testset "Half" begin
    rng = MersenneTwister(42)
    
    @testset "Basic properties" begin
        μ = Half(StdNormal())
        
        # Test show method
        @test sprint(show, μ) == "Half(StdNormal())"
        
        # Test unhalf
        @test unhalf(μ) === StdNormal()
        
        # Test basemeasure
        @test basemeasure(μ) isa WeightedMeasure
        @test _logweight(basemeasure(μ)) ≈ log(2)
    end

    @testset "Sampling and density" begin
        μ = Half(StdNormal())
        n_samples = 1000
        samples = [rand(rng, μ) for _ in 1:n_samples]
        
        # All samples should be non-negative
        @test all(x -> x ≥ 0, samples)
        
        # Test density at specific points
        x = 1.0
        expected_log_density = -0.5 * (x^2 + log2π) - loghalf
        @test logdensityof(μ, x) ≈ expected_log_density
        
        # Test density at negative points
        @test logdensityof(μ, -1.0) == -Inf
        
        # Test density at zero
        @test isfinite(logdensityof(μ, 0.0))
    end

    @testset "Transport" begin
        μ = Half(StdNormal())
        
        # Test transport to/from uniform
        u = 0.7  # arbitrary point in (0,1)
        x = transport_to(μ, StdUniform(), u)
        @test x ≥ 0
        @test transport_to(StdUniform(), μ, x) ≈ u
        
        # Test edge cases
        @test transport_to(μ, StdUniform(), 0.0) == 0.0
        @test transport_to(μ, StdUniform(), 1.0) > 0
    end

    @testset "SMF (Standardized Measure Function)" begin
        μ = Half(StdUniform())
        
        # Test SMF properties
        @test smf(μ, -1.0) == -1.0  # Below support
        @test smf(μ, 0.0) == -1.0    # At lower bound
        @test smf(μ, 0.5) == 0.0     # Midpoint
        @test smf(μ, 1.0) == 1.0     # At upper bound
        
        # Test inverse SMF
        for p in [0.0, 0.25, 0.5, 0.75, 1.0]
            x = invsmf(μ, p)
            @test smf(μ, x) ≈ p
            @test 0 ≤ x ≤ 1
        end
    end
end