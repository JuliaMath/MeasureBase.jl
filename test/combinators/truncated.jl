using Test
using StableRNGs: StableRNG

using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential, StdLogistic, TruncatedMeasure, NormalizedMeasure
using MeasureBase: normalize, massof, smf, invsmf, logdensities, transport_to, mspace_ndims, getdof, fast_dof
using MeasureBase.Interface: test_transport, test_smf
using ChangesOfVariables: with_logabsdet_jacobian
using InverseFunctions: inverse
using IntervalSets
import JLArrays
using JLArrays: JLArray

# Composite Simpson rule:
function _simpson(f, a, b, n = 4000)
    h = (b - a) / n
    h / 3 * (f(a) + f(b) + 4 * sum(i -> f(a + (2i - 1) * h), 1:(n ÷ 2)) + 2 * sum(i -> f(a + 2i * h), 1:(n ÷ 2 - 1)))
end

@testset "truncate" begin
    # Truncations with finite integration ranges (tails beyond them are
    # negligible):
    cases = [
        (StdNormal(), -1.0 .. 2.0, (-1.0, 2.0)),
        (StdNormal(), 0.5 .. Inf, (0.5, 40.0)),
        (StdLogistic(), -Inf .. 0.3, (-60.0, 0.3)),
        (StdExponential(), 0.5 .. 3, (0.5, 3.0)),
        (2.5 * StdNormal(), -0.3 .. 1.2, (-0.3, 1.2)),
        (Half(StdNormal()), 0.5 .. 1.0, (0.5, 1.0)),
        (Lebesgue(), 0 .. 2, (0.0, 2.0)),
    ]

    @testset "masses and densities of truncate($μ, $a_b)" for (μ, a_b, (a, b)) in cases
        t = @inferred truncate(μ, a_b)
        @test t isa TruncatedMeasure
        @test mspace_ndims(t) == 0 && getdof(t) == 1 && fast_dof(t) == 1
        density(x) = exp(logdensityof(μ, x))
        @test massof(t) ≈ _simpson(density, a, b) rtol = 1e-10
        @test massof(t) ≈ massof(μ, a_b)

        xs = range(a, b, length = 7)
        @test logdensityof.(Ref(t), xs) == logdensityof.(Ref(μ), xs)
        @test logdensities(t, collect(xs)) ≈ logdensityof.(Ref(μ), xs)
        lo, hi = float.(endpoints(a_b))
        for x in [filter(isfinite, [prevfloat(lo), nextfloat(hi)]); NaN]
            @test logdensityof(t, x) == -Inf
            @test !insupportof(t, x)
        end

        n = normalize(t)
        @test n isa NormalizedMeasure && massof(n) == 1
        @test _simpson(x -> exp(logdensityof(n, x)), a, b) ≈ 1 rtol = 1e-10
    end

    @testset "open and closed intervals" begin
        @test truncate(StdNormal(), OpenInterval(-1.0, 2.0)) == truncate(StdNormal(), -1.0 .. 2.0)
        @test truncate(StdNormal(), Interval{:closed,:open}(-1.0, 2.0)) == truncate(StdNormal(), -1.0 .. 2.0)
    end

    @testset "invalid truncations" begin
        @test_throws ArgumentError truncate(StdNormal()^2, 0 .. 1)
        @test_throws ArgumentError truncate(MeasureBase.Dirac(0.5), 0 .. 1)
        @test_throws ArgumentError truncate(Counting(Int), 0 .. 1)
        @test_throws ArgumentError truncate(restrict(x -> x > 0, StdNormal()), 0 .. 1)
        @test_throws ArgumentError truncate(StdNormal(), 2.0 .. 1.0)
        @test_throws ArgumentError normalize(truncate(StdNormal(), 1.0 .. 1.0))
        t = truncate(StdNormal(), -1.0 .. 2.0)
        @test_throws ArgumentError logdensityof(t, [0.5])
        @test_throws ArgumentError transport_to(StdUniform(), t)([0.5])
    end

    @testset "transports of truncate($μ, $a_b)" for (μ, a_b, (a, b)) in cases
        t = truncate(μ, a_b)
        n = normalize(t)
        lo, hi = float.(endpoints(a_b))
        f = transport_to(StdUniform(), t)
        # The derivative of the conditional Stieltjes measure function is the
        # normalized density:
        for u in 0.1:0.1:0.9
            x = inverse(f)(u)
            @test a < x < b
            @test f(x) ≈ u
            h = 1e-5 * max(1, abs(x))
            @test (f(x + h) - f(x - h)) / 2h ≈ exp(logdensityof(n, x)) rtol = 1e-6
            @test last(with_logabsdet_jacobian(f, x)) ≈ logdensityof(n, x)
        end
        @test f(lo) == 0
        @test f(hi) ≈ 1
        @test all(isnan, f.([filter(isfinite, [prevfloat(lo), nextfloat(hi)]); NaN]))
        # Finite inputs give finite variates within the truncation:
        for u in (0.0, 1.0, nextfloat(0.0), prevfloat(1.0))
            x = inverse(f)(u)
            @test isfinite(x) && lo <= x <= hi
        end
        @test isnan(inverse(f)(-0.1)) && isnan(inverse(f)(1.1))

        test_transport(StdNormal(), t)
        test_transport(t, StdLogistic())
        test_smf(n)
        isfinite(lo) && @test invsmf(n, 0.0) ≈ lo
        isfinite(hi) && @test invsmf(n, 1.0) ≈ hi

        X = rand(StableRNG(789990641), t^(3, 4))
        @test f.(X) ≈ map(f, X)
        @test inverse(f).(f.(X)) ≈ X
    end

    @testset "random variates of truncate($μ, $a_b)" for (μ, a_b, (a, b)) in cases
        t = truncate(μ, a_b)
        n = normalize(t)
        m1 = _simpson(x -> x * exp(logdensityof(n, x)), a, b)
        m2 = _simpson(x -> (x - m1)^2 * exp(logdensityof(n, x)), a, b)
        N = 10^5
        X = rand(StableRNG(789990641), t^N)
        @test all(x -> a <= x <= b, X)
        @test sum(X) / N ≈ m1 atol = 5 * sqrt(m2 / N)
        @test sum(x -> (x - m1)^2, X) / N ≈ m2 rtol = 0.03
        # Empirical cdf:
        for x in range(a, b, length = 5)[2:(end - 1)]
            @test count(<=(x), X) / N ≈ smf(n, x) atol = 0.01
        end
        x = rand(StableRNG(789990641), t)
        @test x isa Float64 && a <= x <= b
    end

    @testset "precision" begin
        t = truncate(StdNormal(), -1.0 .. 2.0)
        @test rand(StableRNG(1), Float32, t) isa Float32
        @test eltype(rand(StableRNG(1), Float32, t^10)) == Float32
        t32 = truncate(StdNormal(), -1.0f0 .. 2.0f0)
        @test massof(t32) isa Float32
        f32 = transport_to(StdUniform(), t32)
        @test f32(0.5f0) isa Float32 && inverse(f32)(0.5f0) isa Float32
        @test f32(0.5f0) ≈ transport_to(StdUniform(), t)(0.5) rtol = 1e-6
    end

    @testset "device arrays" begin
        JLArrays.allowscalar(false)
        t = truncate(StdNormal(), -1.0 .. 2.0)
        x = [-2.0, 0.0, 0.5, 1.9, 3.0]
        @test all(Array(logdensities(t, JLArray(x))) .≈ logdensities(t, x))
        f = transport_to(StdUniform(), t)
        @test isequal(Array(f.(JLArray(x))), f.(x))
        @test Array(inverse(f).(JLArray([0.0, 0.3, 1.0]))) ≈ inverse(f).([0.0, 0.3, 1.0])
    end

    @testset "Half and normalized truncations" begin
        t = normalize(truncate(StdNormal(), 0 .. Inf))
        h = Half(StdNormal())
        for x in (0.0, 0.3, 2.5)
            @test logdensityof(t, x) ≈ logdensityof(h, x)
            @test transport_to(StdUniform(), t)(x) ≈ transport_to(StdUniform(), h)(x)
        end
        @test logdensityof(t, -0.3) == logdensityof(h, -0.3) == -Inf
        # Truncations of weighted measures keep the weight, their transports
        # don't depend on it:
        tw = truncate(2.5 * StdNormal(), -1.0 .. 2.0)
        t1 = truncate(StdNormal(), -1.0 .. 2.0)
        @test massof(tw) ≈ 2.5 * massof(t1)
        @test transport_to(StdNormal(), tw)(0.3) ≈ transport_to(StdNormal(), t1)(0.3)
        @test logdensityof(normalize(tw), 0.3) ≈ logdensityof(normalize(t1), 0.3)
        # Nested truncations:
        tt = truncate(t1, 0.0 .. 3.0)
        @test massof(tt) ≈ massof(truncate(StdNormal(), 0.0 .. 2.0))
        @test transport_to(StdUniform(), tt)(1.2) ≈ transport_to(StdUniform(), truncate(StdNormal(), 0.0 .. 2.0))(1.2)
    end
end
