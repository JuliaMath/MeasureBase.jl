# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: smf, invsmf, massof, normalize, transport_to, logdensities
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
