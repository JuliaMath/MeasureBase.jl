# test/domains.jl
using Test
using MeasureBase
using Static: static
using Random: MersenneTwister

@testset "Domains" begin
    @testset "BoundedInts" begin
        bounded = ℤ[1:5]
        @test 3 ∈ bounded
        @test -1 ∉ bounded
        @test 6 ∉ bounded
        @test 1.5 ∉ bounded
        @test minimum(bounded) == 1
        @test maximum(bounded) == 5
        @test testvalue(bounded) == 0
        
        # Test show method
        @test sprint(show, bounded) == "ℤ[1:5]"
    end

    @testset "ZeroSet" begin
        # Simple quadratic function and its gradient
        f(x) = sum(x.^2)
        ∇f(x) = 2x
        zs = ZeroSet(f, ∇f)
        
        # Test points
        @test zeros(3) ∈ zs
        @test [1e-8, -1e-8, 1e-8] ∈ zs
        @test [0.1, 0.1, 0.1] ∉ zs
        
        # Test with different floating point types
        @test zeros(Float32, 2) ∈ zs
        @test zeros(Float64, 2) ∈ zs
    end

    @testset "IntegerNumbers" begin
        @test minimum(ℤ) == static(-Inf)
        @test maximum(ℤ) == static(Inf)
        
        # Test membership
        @test 42 ∈ ℤ
        @test -42 ∈ ℤ
        @test 3.14 ∉ ℤ
        @test 2.0 ∈ ℤ  # Integer-valued floats should be in ℤ
    end
end