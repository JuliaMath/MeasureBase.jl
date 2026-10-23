# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: smf, invsmf, massof, normalize, transport_to, logdensities
using MeasureBase: csmf, invcsmf, logsmf, logcsmf, invlogsmf, invlogcsmf
using Distributions
using IntervalSets
using InverseFunctions: inverse
using StableRNGs: StableRNG
import ForwardDiff

@testset "Stieltjes measure functions of distributions" begin
    for d in [Normal(0.3, 1.7), Gamma(2.0, 1.5), Beta(2.0, 3.0), Weibull(0.7, 1.3)]
        m = asmeasure(d)
        for x in quantile.(d, [0.01, 0.3, 0.7, 0.99])
            @test smf(m, x) ≈ cdf(d, x)
            @test invsmf(m, cdf(d, x)) ≈ x
        end
        @test massof(m, quantile(d, 0.2) .. quantile(d, 0.7)) ≈ 0.5
        @test isnan(smf(m, NaN)) && isnan(invsmf(m, NaN)) && isnan(invsmf(m, 1.5))
    end
    @test smf(asmeasure(Poisson(2.0)), 2.5) ≈ cdf(Poisson(2.0), 2)
    @test smf(asmeasure(Normal()), Inf) == 1 && smf(asmeasure(Normal()), -Inf) == 0
    @test smf(3 * asmeasure(Normal()), 0.0) ≈ 1.5
end

@testset "complementary and log Stieltjes measure functions of distributions" begin
    for d in [Normal(0.3, 1.7), Gamma(2.0, 1.5), Beta(2.0, 3.0), Weibull(0.7, 1.3)]
        m = asmeasure(d)
        # Both tails, beyond the resolution of the cdf resp. ccdf near one:
        for p in [1e-30, 1e-19, 0.01, 0.3, 0.7]
            for x in (quantile(d, p), cquantile(d, p))
                @test csmf(m, x) ≈ ccdf(d, x)
                @test logcsmf(m, x) ≈ logccdf(d, x)
                @test logsmf(m, x) ≈ logcdf(d, x)
            end
            @test invcsmf(m, p) ≈ cquantile(d, p)
            @test invlogcsmf(m, log(p)) ≈ invlogccdf(d, log(p))
            @test invlogsmf(m, log(p)) ≈ invlogcdf(d, log(p))
        end
        @test isnan(csmf(m, NaN)) && isnan(logcsmf(m, NaN)) && isnan(invcsmf(m, 1.5)) && isnan(invlogcsmf(m, 0.5))
    end
    m = asmeasure(Normal())
    @test csmf(m, Inf) == 0 && csmf(m, -Inf) == 1
    @test logcsmf(m, Inf) == -Inf && logcsmf(m, -Inf) == 0 && logsmf(m, -Inf) == -Inf && logsmf(m, Inf) == 0
    @test csmf(asmeasure(Poisson(2.0)), 2.5) ≈ ccdf(Poisson(2.0), 2)
    @test csmf(3 * m, 9.0) ≈ 3 * ccdf(Normal(), 9.0)
end

@testset "truncated distributions" begin
    for (d, a, b) in [
        (Normal(1.0, 2.0), -0.5, 2.5),
        (Normal(1.0, 2.0), 0.0, Inf),
        (Gamma(2.0, 1.5), 1.0, 4.0),
        (Exponential(2.0), -Inf, 3.0),
        (Beta(2.0, 3.0), 0.2, 0.9),
        (Logistic(0.5, 2.0), -Inf, 0.0),
    ]
        t = truncate(asmeasure(d), a .. b)
        n = normalize(t)
        dref = truncated(d, a, b)
        @test massof(t) ≈ cdf(d, b) - cdf(d, a)
        xs = [quantile.(dref, [0.01, 0.2, 0.5, 0.8, 0.99]); a - 1; b + 1]
        @test all(logdensityof.(Ref(n), xs) .≈ logpdf.(dref, xs))
        @test all(logdensities(n, xs) .≈ logpdf.(dref, xs))
        f = transport_to(StdUniform(), t)
        @test f.(xs[1:5]) ≈ [0.01, 0.2, 0.5, 0.8, 0.99]
        @test inverse(f).([0.01, 0.2, 0.5, 0.8, 0.99]) ≈ xs[1:5]
        @test all(isnan, f.(xs[6:7]))
        X = rand(StableRNG(789990641), n^(10^5))
        @test minimum(dref) <= minimum(X) && maximum(X) <= maximum(dref)
        # Moments by midpoint quantile quadrature:
        Q = quantile.(dref, ((1:10^5) .- 0.5) ./ 10^5)
        @test sum(X) / length(X) ≈ mean(Q) atol = 5 * std(Q) / sqrt(length(X))
    end

    # Derivatives with respect to the parameters, also for half-infinite
    # intervals:
    ℓ(μ, x) = logdensityof(normalize(truncate(asmeasure(Normal(μ, 1.5)), 0.0 .. Inf)), x)
    ℓ_ref(μ, x) = logpdf(Normal(μ, 1.5), x) - logccdf(Normal(μ, 1.5), 0.0)
    for x in (0.2, 1.7)
        @test ℓ(0.4, x) ≈ ℓ_ref(0.4, x)
        @test ForwardDiff.derivative(μ -> ℓ(μ, x), 0.4) ≈ ForwardDiff.derivative(μ -> ℓ_ref(μ, x), 0.4)
    end
    u(μ) = transport_to(StdUniform(), truncate(asmeasure(Normal(μ, 1.5)), -1.0 .. Inf))(0.7)
    @test ForwardDiff.derivative(u, 0.4) ≈ (u(0.4 + 1e-6) - u(0.4 - 1e-6)) / 2e-6 rtol = 1e-6

    @test rand(StableRNG(1), Float32, truncate(asmeasure(Normal(1.0, 2.0)), 0.0 .. Inf)) isa Float32
    @test_throws ArgumentError truncate(asmeasure(Poisson(2.0)), 0 .. 3)
    @test_throws ArgumentError truncate(asmeasure(MvNormal([1.0 0.0; 0.0 1.0])), 0 .. 3)
end

@testset "truncations to tails" begin
    # Masses, densities and transports of truncations to far tails, against
    # the ccdf in the upper and the cdf in the lower tail:
    for (μ, d) in [(StdNormal(), Normal()), (asmeasure(Normal(0.3, 1.7)), Normal(0.3, 1.7)), (asmeasure(Gamma(2.0, 1.5)), Gamma(2.0, 1.5))]
        for (a, b, upper) in [
            (cquantile(d, 1e-19), Inf, true),
            (cquantile(d, 1e-19), cquantile(d, 1e-22), true),
            (minimum(d), quantile(d, 1e-19), false),
        ]
            t = truncate(μ, a .. b)
            mass_ref = upper ? ccdf(d, a) - ccdf(d, b) : cdf(d, b) - cdf(d, a)
            @test massof(t) ≈ mass_ref rtol = 1e-10
            @test massof(μ, a .. b) ≈ mass_ref rtol = 1e-10
            n = normalize(t)
            f = transport_to(StdUniform(), t)
            for u in (1e-6, 0.1, 0.5, 0.9)
                x = upper ? cquantile(d, ccdf(d, a) - u * mass_ref) : quantile(d, cdf(d, a) + u * mass_ref)
                @test f(x) ≈ u rtol = 1e-6
                @test inverse(f)(u) ≈ x rtol = 1e-12
                @test logdensityof(n, x) ≈ logpdf(d, x) - log(mass_ref)
                @test smf(t, x) + csmf(t, x) ≈ massof(t)
            end
            @test all(isfinite, inverse(f).([0.0, 1.0]))
        end
    end

    t32 = truncate(StdNormal(), 9.0f0 .. Inf32)
    @test massof(t32) isa Float32
    @test massof(t32) ≈ ccdf(Normal(), 9.0) rtol = 1e-5
    f32 = transport_to(StdUniform(), t32)
    @test f32(9.5f0) isa Float32
    @test f32(9.5f0) ≈ 1 - ccdf(Normal(), 9.5) / ccdf(Normal(), 9.0) rtol = 1e-5
    @test inverse(f32)(0.5f0) isa Float32
    @test inverse(f32)(0.5f0) ≈ cquantile(Normal(), ccdf(Normal(), 9.0) / 2) rtol = 1e-5
    @test isfinite(inverse(f32)(1.0f0))
end
