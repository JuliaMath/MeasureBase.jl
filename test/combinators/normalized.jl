using Test
using Random: MersenneTwister

using MeasureBase
using MeasureBase: StdNormal, StdExponential, StdUniform, NormalizedMeasure, BoundedInts
using MeasureBase: normalize, massof, isnormalized, UnknownMass, mbind, weightedmeasure, logdensities
using MeasureBase.Interface: test_transport
using ChangesOfVariables: with_logabsdet_jacobian
using InverseFunctions: inverse
using AffineMaps: Add, Mul
using IntervalSets
using Static: static

# Midpoint rule on a uniform grid, exact to O(h²) for piecewise smooth
# integrands with jumps on grid points:
_midpoint_integral(f, a, b, n) = (h = (b - a) / n; h * sum(i -> f(a + (i - 0.5) * h), 1:n))

# Twice the standard exponential measure, with a transport, with and
# without a declared mass:
struct _TwiceStdExp{DeclaredMass} <: MeasureBase.AbstractMeasure end
MeasureBase.basemeasure(::_TwiceStdExp) = StdExponential()
MeasureBase.logdensity_def(::_TwiceStdExp, x) = log(2.0)
MeasureBase.insupportof(::_TwiceStdExp, x) = x >= 0
MeasureBase.mspace_flatsize(::Type{<:_TwiceStdExp}) = ()
MeasureBase.preferred_stdmeasure(::Type{<:_TwiceStdExp}) = StdExponential
MeasureBase.transport_to_std(::Type{StdExponential}, ::_TwiceStdExp, x) = x
MeasureBase.transport_from_std(::Type{StdExponential}, ::_TwiceStdExp, z) = z
MeasureBase.massof(::_TwiceStdExp{true}) = 2.0

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

    @testset "transport Jacobians" begin
        # Array products are normalized via their mass, the Jacobian of the
        # transport to the uniform distribution is the normal density:
        P = MeasureBase.ProductMeasure([weightedmeasure(log(i), StdNormal()) for i in 1:3])
        x = [0.3, -1.2, 2.0]
        @test last(with_logabsdet_jacobian(transport_to(StdUniform()^3, P), x)) ≈ sum(@. -x^2 / 2 - log(2π) / 2)
        test_transport(StdUniform()^3, normalize(P))
        test_transport(normalize(P), StdNormal()^3)
        # Powers of products of weighted measures, z ~ N(0, 1) to
        # (normal, exponential) pairs:
        Q = productmeasure((2 * StdNormal(), StdExponential()))^2
        g = transport_to(Q, StdNormal()^4)
        z = [0.3, -0.5, 1.1, 0.8]
        y = g(z)
        ℓ_trg = sum(t -> -t[1]^2 / 2 - log(2π) / 2 - t[2], y)
        @test last(with_logabsdet_jacobian(g, z)) ≈ sum(@. -z^2 / 2 - log(2π) / 2) - ℓ_trg
        # Pushforwards of binds of unknown mass are normalized structurally,
        # kernels weighted or not. Variates [a, y] with a ~ U(0, 1),
        # y ~ N(0, 1) resp. N(a, 1), pushed forward by Mul(2) and transported
        # to [Φ⁻¹(a), y] resp. [Φ⁻¹(a), y - a]:
        wb = mbind(a -> weightedmeasure(only(a), StdNormal()^1), StdUniform()^1, vcat)
        nb = mbind(a -> pushfwd(Add(only(a)), StdNormal()^1), StdUniform()^1, vcat)
        for β in (wb, nb)
            h = transport_to(StdNormal()^2, pushfwd(Mul(2.0), β))
            @test h([0.8, 0.4]) ≈ transport_to(StdNormal()^2, β)([0.4, 0.2])
            @test last(with_logabsdet_jacobian(h, [0.8, 0.4])) ≈ -2 * log(2) - logdensityof(StdNormal(), h([0.8, 0.4])[1])
        end
        # u = 1 - exp(-x) has the log-Jacobian -x if the mass is known,
        # without a known mass the Jacobian has no normalized density:
        t = transport_to(StdUniform(), _TwiceStdExp{true}())
        @test last(with_logabsdet_jacobian(t, 1.2)) ≈ -1.2
        @test last(with_logabsdet_jacobian(inverse(t), t(1.2))) ≈ 1.2
        u = transport_to(StdUniform(), _TwiceStdExp{false}())
        @test u(1.2) ≈ t(1.2)
        @test_throws ArgumentError with_logabsdet_jacobian(u, 1.2)
        @test_throws ArgumentError with_logabsdet_jacobian(inverse(u), u(1.2))
    end
end
