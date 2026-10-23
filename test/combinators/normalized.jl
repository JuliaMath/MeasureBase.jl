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

# Midpoint rule on a uniform grid, exact to O(h²) for piecewise smooth
# integrands with jumps on grid points:
_midpoint_integral(f, a, b, n) = (h = (b - a) / n; h * sum(i -> f(a + (i - 0.5) * h), 1:n))

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

    s = superpose(0.5 * StdNormal(), 2 * StdExponential())
    b = mbind(a -> pushfwd(Add(a), StdNormal()), StdNormal(), tuple)

    @testset "construction" begin
        @test @inferred(normalize(StdNormal())) === StdNormal()
        @test normalize(Dirac(1.5)) === Dirac(1.5)
        @test @inferred(normalize(3 * StdNormal())) === StdNormal()
        @test @inferred(normalize(s)) isa NormalizedMeasure{typeof(s),Float64}
        @test @inferred(normalize(b)) isa NormalizedMeasure{typeof(b),UnknownMass}
        for μ in (StdNormal(), Lebesgue(0 .. 2), s, b, Counting(BoundedInts(0, 10)))
            @test normalize(normalize(μ)) === normalize(μ)
            @test normalize(weightedmeasure(0.7, μ)) == normalize(μ)
            @test massof(normalize(μ)) == 1 && isnormalized(normalize(μ))
        end
        @test_throws ArgumentError normalize(Lebesgue())
        @test_throws ArgumentError normalize(Lebesgue(2 .. 1))
        @test_throws ArgumentError normalize(Counting())
    end

    @testset "densities" begin
        @test logdensityof(normalize(3 * StdNormal()), 0.7) == logdensityof(StdNormal(), 0.7)
        u = normalize(Lebesgue(0 .. 2))
        @test @inferred(logdensityof(u, 0.5)) ≈ log(1 / 2)
        @test logdensityof(u, 2.5) == -Inf
        @test logdensity_rel(u, Lebesgue(0 .. 2), 1.0) ≈ -log(2)
        ν = normalize(s)
        @test @inferred(logdensityof(ν, 0.3)) ≈ log((0.5 * exp(logdensityof(StdNormal(), 0.3)) + 2 * exp(-0.3)) / 2.5)
        @test _midpoint_integral(x -> exp(logdensityof(ν, x)), -12.0, 40.0, 52_000) ≈ 1 atol = 1e-6
        c = normalize(Counting(BoundedInts(0, 10)))
        @test sum(k -> exp(logdensityof(c, k)), -2:12) ≈ 1
        d = normalize(superpose(2 * Dirac(1.0), 3 * Dirac(2.0)))
        @test exp(logdensityof(d, 1.0)) ≈ 2 / 5 && exp(logdensityof(d, 2.0)) ≈ 3 / 5
    end

    @testset "batched densities" begin
        X = [-1.0, 0.3, 2.0]
        @test logdensities(normalize(s), X) ≈ logdensities(s, X) .- log(2.5)
        @test logdensities(normalize(2.0 * Lebesgue(0 .. 2)), [0.5, 3.0]) == [-log(2), -Inf]
        P = superpose(StdNormal()^2, 2 * StdNormal()^2)
        Y = randn(2, 4)
        @test logdensities(normalize(P), Y) ≈ logdensities(P, Y) .- log(3)
        @test logdensities(normalize(Lebesgue(0 .. 2))^2, [0.5 1.0; 1.5 2.5]) ≈ [-2 * log(2), -Inf]
    end

    @testset "unknown mass" begin
        ν = normalize(b)
        x = rand(ν)
        @test x isa Tuple{Float64,Float64}
        @test rand(MersenneTwister(7), ν) == rand(MersenneTwister(7), b)
        @test insupportof(ν, x) == insupportof(b, x)
        @test_throws ArgumentError logdensityof(ν, x)
        @test_throws ArgumentError logdensities(ν, [x, x])
        # The transport of the bind, z = (a, y - a), has unit Jacobian:
        f = transport_to(StdNormal()^2, ν)
        @test f((0.3, 1.0)) ≈ [0.3, 0.7]
        @test last(with_logabsdet_jacobian(f, (0.3, 1.0))) ≈ 0 atol = 1e-12
    end
end
