using Test
using Random: MersenneTwister

using MeasureBase
using MeasureBase: StdNormal, StdExponential, StdUniform, NormalizedMeasure, BoundedInts
using MeasureBase: normalize, massof, isnormalized, UnknownMass, mbind, weightedmeasure, logdensities
using MeasureBase.Interface: test_transport
using ChangesOfVariables: with_logabsdet_jacobian
using AffineMaps: Add, Mul
using IntervalSets
using Static: static

@testset "normalize" begin
    @testset "known masses" begin
        @test @inferred(massof(Lebesgue(0 .. 2))) == 2
        @test massof(Lebesgue(1.0 .. 3.5)) == 2.5
        @test @inferred(massof(Counting(BoundedInts(0, 10)))) == 11
        @test massof(Counting(BoundedInts(3, 2))) == 0
        @test @inferred(massof(Counting(BoundedInts(static(1), static(4))))) === static(4)
        @test massof(Counting(BoundedInts(static(0), static(Inf)))) == Inf
        @test massof(Counting(1:4)) == 4 && massof(Counting(Set([1, 3]))) == 2
        @test massof(Counting()) == Inf
        @test @inferred(massof(3 * StdNormal())) ≈ 3
        @test @inferred(massof(productmeasure((2 * StdNormal(), Lebesgue(0 .. 3))))) ≈ 6
        @test @inferred(massof(Lebesgue(0 .. 2)^3)) == 8
        @test @inferred(massof(superpose(StdNormal(), 2 * Lebesgue(0 .. 2)))) ≈ 5
        @test @inferred(massof(StdNormal()^3)) === static(1.0)
        @test @inferred(massof(productmeasure((StdNormal(), restrict(x -> x > 0, StdNormal()))))) isa UnknownMass
    end
end
