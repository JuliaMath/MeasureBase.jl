# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test
using MeasureBase
using MeasureBase: smf, csmf, logsmf, logcsmf, invsmf, invcsmf, invlogsmf, invlogcsmf
using Distributions
import ForwardDiff, Enzyme

@testset "affine Stieltjes measure functions" begin
    fd(f, x) = (f(x + 1e-5) - f(x - 1e-5)) / 2e-5

    for σ in (3.0, -3.0)
        μ, α, z = 1.0, 2.5, 0.4
        x = μ + σ * z
        m(a) = asmeasure(μ + σ * Beta(a, 1.0))
        p = σ > 0 ? z^α : 1 - z^α
        dp = (σ > 0 ? 1 : -1) * z^α * log(z)
        for (f, value, deriv) in (
            (smf, p, dp),
            (csmf, 1 - p, -dp),
            (logsmf, log(p), dp / p),
            (logcsmf, log1p(-p), -dp / (1 - p)),
        )
            @test f(m(α), x) ≈ value
            @test ForwardDiff.derivative(a -> f(m(a), x), α) ≈ deriv rtol = 1e-10
        end

        for (f, fc) in ((smf, csmf), (logsmf, logcsmf))
            expected = f === smf ? (0.0, 1.0) : (-Inf, 0.0)
            @test (f(m(α), -Inf), f(m(α), Inf)) == expected
            @test (fc(m(α), Inf), fc(m(α), -Inf)) == expected
            @test ForwardDiff.derivative(s -> f(asmeasure(μ + s * Beta(α, 1.0)), Inf), σ) ==
                  0
            @test ForwardDiff.derivative(
                l -> fc(asmeasure(l + σ * Beta(α, 1.0)), -Inf),
                μ,
            ) == 0
        end

        for p in (0.0, 0.3, 1.0),
            (f, complement, logarithmic) in (
                (invsmf, false, false),
                (invcsmf, true, false),
                (invlogsmf, false, true),
                (invlogcsmf, true, true),
            )

            u = xor(complement, σ < 0) ? 1 - p : p
            arg = logarithmic ? log(p) : p
            ref = μ + σ * u^(1 / α)
            deriv = 0 < u < 1 ? -σ * u^(1 / α) * log(u) / α^2 : 0.0
            @test f(m(α), arg) ≈ ref
            @test ForwardDiff.derivative(a -> f(m(a), arg), α) ≈ deriv rtol = 1e-9 atol =
                1e-14
        end

        q(a) = invsmf(m(a), 0.3)
        @test ForwardDiff.derivative(p -> invsmf(m(α), p), 0.3) ≈
              1 / pdf(μ + σ * Beta(α, 1.0), q(α))
        @test ForwardDiff.derivative(a -> smf(m(a), q(a)), α) ≈ 0 atol = 1e-14
        @test ForwardDiff.derivative(l -> smf(asmeasure(l + σ * Beta(α, 1.0)), x), μ) ≈
              -pdf(μ + σ * Beta(α, 1.0), x)
        @test ForwardDiff.derivative(s -> smf(asmeasure(μ + s * Beta(α, 1.0)), x), σ) ≈
              -z * pdf(μ + σ * Beta(α, 1.0), x)
        a = min(μ, μ + σ) + 0.2
        interval_q(α) = invsmf(m(α), a, 0.4)
        @test smf(m(α), a, interval_q(α)) ≈ 0.4
        @test ForwardDiff.derivative(interval_q, α) ≈ fd(interval_q, α) rtol = 1e-8
    end

    f(p) =
        smf(asmeasure(LocationScale(122.0, 6.0, Beta(exp(p[1]), exp(p[2])))), 124.0, 126.0)
    p = [log(12.0), log(12.0)]
    reference = [fd(a -> f([a, p[2]]), p[1]), fd(b -> f([p[1], b]), p[2])]
    @test ForwardDiff.gradient(f, p) ≈ reference rtol = 1e-8

    for f in (invsmf, invcsmf)
        g(α) = f(asmeasure(1.0 + (-2.0) * Gamma(α, 1.5)), 0.3)
        @test ForwardDiff.derivative(g, 2.5) ≈ fd(g, 2.5) rtol = 1e-8
    end
    for σ in (2.0, -2.0)
        m = asmeasure(1.0 + σ * Gamma(2.0, 1.5))
        f, g = σ > 0 ? (csmf, invcsmf) : (smf, invsmf)
        @test f(m, g(m, 1e-30)) ≈ 1e-30 rtol = 1e-12
    end
    @test invcsmf(asmeasure(1.0 + 2.0 * Exponential()), 0.0) ≈ 1 - 2log(floatmin(Float64))

    m32(a) = asmeasure(1.0f0 + 3.0f0 * Beta(a, 1.0f0))
    @test smf(m32(2.0f0), 2.5f0) ≈ 0.25f0
    @test ForwardDiff.derivative(a -> smf(m32(a), 2.5f0), 2.0f0) ≈ 0.25f0 * log(0.5f0) rtol =
        1e-5
    @test invsmf(m32(2.0f0), 0.25f0) ≈ 2.5f0
    @test ForwardDiff.derivative(p -> invsmf(m32(2.0f0), p), 0.25f0) ≈ 3.0f0 rtol = 1e-5

    # The small input tail must survive even when the beta variate is above 1/2.
    p = 1.0f-6
    x = exp(log(Float64(p)) / 50)
    y = -expm1(log(Float64(p)) / 50)
    for (law, expected) in (
        (Beta(50.0f0, 1.0f0), x),
        (Distributions.BetaPrime(50.0f0, 1.0f0), x / y),
        (Distributions.FDist(100.0f0, 2.0f0), x / (50y)),
    )
        @test invsmf(asmeasure(1.0f0 + 2.0f0 * law), p) ≈ 1 + 2expected rtol = 8eps(Float32)
    end

    # The beta coordinate rounds to one here; its complement still gives a finite ratio.
    point = [0.02f0, 0.3f0]
    b, p = Float64.(point)
    w = -log1p(-p) / b
    expected = [exp(w) * log1p(-p) / b^2, exp(w) / (b * (1-p))]
    quantile_ratio(v) = invsmf(asmeasure(Distributions.BetaPrime(1.0f0, v[1])), v[2])
    @test quantile_ratio(point) ≈ expm1(w) rtol = 32eps(Float32)
    @test ForwardDiff.gradient(quantile_ratio, point) ≈ expected rtol = 64eps(Float32)
    @test only(Enzyme.gradient(Enzyme.Reverse, quantile_ratio, point)) ≈ expected rtol =
        64eps(Float32)
end
