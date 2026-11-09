# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: smf, invsmf, massof, normalize, transport_to, logdensities, StdLogistic
using MeasureBase: csmf, invcsmf, logsmf, logcsmf, invlogsmf, invlogcsmf
using Distributions
using IntervalSets
using InverseFunctions: inverse
using StableRNGs: StableRNG
using SpecialFunctions: erf, erfc
import ForwardDiff, Zygote, Mooncake

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

@testset "Cauchy Stieltjes measure functions in the tails" begin
    c = asmeasure(Cauchy())
    for x in (1e15, 1e17, 1e300)
        @test csmf(c, x) ≈ 1 / (π * x) && smf(c, -x) ≈ 1 / (π * x)
        @test logcsmf(c, x) ≈ -log(π * x) && logsmf(c, -x) ≈ -log(π * x)
        @test invcsmf(c, csmf(c, x)) ≈ x && invsmf(c, smf(c, -x)) ≈ -x
        @test invlogcsmf(c, logcsmf(c, x)) ≈ x && invlogsmf(c, logsmf(c, -x)) ≈ -x
    end
    @test logsmf(c, 1e17) ≈ -1 / (π * 1e17) && logcsmf(c, -1e17) ≈ -1 / (π * 1e17)
    @test smf(c, 0.3) ≈ cdf(Cauchy(), 0.3) && invsmf(c, 0.3) ≈ quantile(Cauchy(), 0.3)
    @test smf(c, Inf) == 1 && csmf(c, Inf) == 0 && logsmf(c, -Inf) == -Inf && isnan(smf(c, NaN))
    @test isnan(invlogsmf(c, NaN)) && isnan(invlogcsmf(c, 0.5)) && isnan(invcsmf(c, 1.5))
    for a in (1e17, -1e17)
        t = truncate(c, a > 0 ? (a .. Inf) : (-Inf .. a))
        @test massof(t) ≈ 1 / (π * 1e17)
        # The tail beyond 2a holds about half of the truncated mass:
        X = rand(StableRNG(789990641), normalize(t)^10^4)
        @test count(x -> abs(x) > 2e17, X) / length(X) ≈ 0.5 atol = 0.03
        @test all(x -> abs(x) >= 1e17, X)
    end
end

@testset "beta tail probabilities at small shapes" begin
    # Beta(a, 1) has the cdf x^a, its upper tail at tiny x isn't small:
    m = asmeasure(Beta(0.01, 1.0))
    x = 1e-20
    p = x^0.01
    @test csmf(m, x) ≈ 1 - p && logcsmf(m, x) ≈ log1p(-p) && logsmf(m, x) ≈ log(p)
    @test transport_to(StdLogistic(), m)(x) ≈ log(p) - log1p(-p)
    q = (2.0^-40)^0.01
    @test transport_to(StdLogistic(), asmeasure(Beta(1.0, 0.01)))(1 - 2.0^-40) ≈ log1p(-q) - log(q)
end

@testset "Stieltjes measure functions of gamma and beta variate families" begin
    # Through the gamma and beta kernels, against Distributions in the bulk
    # and by round trips in the tails:
    for d in [Chisq(2.3), Erlang(3, 1.5), InverseGamma(2.5, 1.3), FDist(3.0, 5.0), BetaPrime(0.4, 0.6)]
        m = asmeasure(d)
        for x in quantile.(d, [0.01, 0.3, 0.7, 0.99])
            @test smf(m, x) ≈ cdf(d, x) && csmf(m, x) ≈ ccdf(d, x)
            @test logsmf(m, x) ≈ logcdf(d, x) && logcsmf(m, x) ≈ logccdf(d, x)
            @test invsmf(m, cdf(d, x)) ≈ x && invcsmf(m, ccdf(d, x)) ≈ x
        end
        for p in [1e-30, 1e-8, 0.3]
            @test smf(m, invsmf(m, p)) ≈ p && csmf(m, invcsmf(m, p)) ≈ p
            @test logsmf(m, invlogsmf(m, log(p))) ≈ log(p) && logcsmf(m, invlogcsmf(m, log(p))) ≈ log(p)
        end
        @test smf(m, 0.0) == smf(m, -1.0) == csmf(m, Inf) == 0 && csmf(m, 0.0) == smf(m, Inf) == 1
        @test logsmf(m, 0.0) == -Inf && logcsmf(m, 0.0) == 0
        @test isnan(smf(m, NaN)) && isnan(logcsmf(m, NaN)) && isnan(invsmf(m, 1.5)) && isnan(invlogcsmf(m, 0.5))
    end
end

@testset "interval masses of distributions" begin
    # μ((a, x]) against references at 256 bits, just above the lower end in
    # the bulk (normal and logistic families via the standard measures) and
    # in the tails, where the Stieltjes measure functions round to one (the
    # gamma and beta families via the complementary forms):
    ref(f, args...) = setprecision(() -> Float64(f(big.(args)...)), BigFloat, 256)
    Φ(z) = erfc(-z / sqrt(big(2))) / 2
    Φ₀(z) = erf(z / sqrt(big(2))) / 2
    σ(z) = 1 / (1 + exp(-z))
    Q2(y) = (1 + y) * exp(-y)
    beta23_ccdf(c) = c^3 * (4 - 3c)
    for (d, a, x, mass) in [
        (Normal(1.0, 2.0), 1.0, 1.0 + 1e-12, (a, x) -> Φ₀((x - 1) / 2) - Φ₀((a - 1) / 2)),
        (Normal(0.0, 2.0), 0.0, 1e-300, (a, x) -> Φ₀(x / 2) - Φ₀(a / 2)),
        (Normal(0.3, 1.7), 15.0, 16.0, (a, x) -> Φ((big(0.3) - a) / big(1.7)) - Φ((big(0.3) - x) / big(1.7))),
        (Normal(0.3, 1.7), -16.0, -15.0, (a, x) -> Φ((x - big(0.3)) / big(1.7)) - Φ((a - big(0.3)) / big(1.7))),
        (Logistic(1.0, 2.0), 1.0, 1.0 + 1e-12, (a, x) -> σ((x - 1) / 2) - σ((a - 1) / 2)),
        (Logistic(0.5, 2.0), 80.0, 90.0, (a, x) -> σ((big(0.5) - a) / 2) - σ((big(0.5) - x) / 2)),
        (Gamma(2.0, 1.5), 60.0, 63.0, (a, x) -> Q2(a / big(1.5)) - Q2(x / big(1.5))),
        (Chisq(4.0), 80.0, 85.0, (a, x) -> Q2(a / 2) - Q2(x / 2)),
        (Erlang(2, 1.5), 60.0, 63.0, (a, x) -> Q2(a / big(1.5)) - Q2(x / big(1.5))),
        (Beta(2.0, 3.0), 0.999, 0.9999, (a, x) -> beta23_ccdf(1 - a) - beta23_ccdf(1 - x)),
        (InverseGamma(2.0, 1.3), 1e3, 1e4, (a, x) -> Q2(big(1.3) / x) - Q2(big(1.3) / a)),
        (FDist(4.0, 6.0), 1e4, 1e5, (a, x) -> beta23_ccdf(6 / (4a + 6)) - beta23_ccdf(6 / (4x + 6))),
        (BetaPrime(2.0, 3.0), 1e3, 1e4, (a, x) -> beta23_ccdf(1 / (1 + a)) - beta23_ccdf(1 / (1 + x))),
    ]
        m = asmeasure(d)
        p = smf(m, a, x)
        @test p ≈ ref(mass, a, x) rtol = 1e-13
        @test smf(d, a, x) == p
        @test invsmf(m, a, p) ≈ x rtol = 1e-13
        @test smf(m, a, a) == 0
    end
end

@testset "shape derivatives of log Stieltjes measure functions" begin
    fd(f, a) = (f(a + 1e-6) - f(a - 1e-6)) / 2e-6
    for D in (a -> Gamma(a, 1.5), a -> Beta(a, 2.0), b -> Beta(1.5, b)), (g, x) in ((logsmf, 0.4), (logcsmf, 0.4), (invlogsmf, -1.2), (invlogcsmf, -1.2), (invlogsmf, -40.0))
        f = a -> g(asmeasure(D(a)), x)
        @test ForwardDiff.derivative(f, 2.5) ≈ fd(f, 2.5) rtol = 1e-6
        @test only(Zygote.gradient(f, 2.5)) ≈ fd(f, 2.5) rtol = 1e-6
    end
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

@testset "reverse-mode derivatives of truncations with infinite ends" begin
    mooncake_derivative(f, x) = Mooncake.value_and_gradient!!(Mooncake.prepare_gradient_cache(f, x), f, x)[2][2]
    ℓ(d, a, b, x) = logdensityof(normalize(truncate(asmeasure(d), a .. b)), x)
    for f in [
        σ -> ℓ(Normal(0.2, σ), 0.0, Inf, 0.5), σ -> ℓ(Normal(0.2, σ), -Inf, 0.0, -0.5),
        σ -> ℓ(LogNormal(0.2, σ), 1.0, Inf, 1.5), σ -> ℓ(Exponential(σ), 1.0, Inf, 1.5),
        σ -> ℓ(Laplace(0.1, σ), 1.0, Inf, 1.5), σ -> ℓ(Logistic(0.1, σ), 1.0, Inf, 1.5),
        σ -> ℓ(Weibull(σ, 1.1), 1.0, Inf, 1.5),
        σ -> transport_to(StdUniform(), truncate(asmeasure(Normal(0.2, σ)), -1.0 .. Inf))(0.7),
    ]
        ∂f_ref = (f(1.3 + 1e-6) - f(1.3 - 1e-6)) / 2e-6
        @test only(Zygote.gradient(f, 1.3)) ≈ ∂f_ref rtol = 1e-6
        @test mooncake_derivative(f, 1.3) ≈ ∂f_ref rtol = 1e-6
    end
    g = σ -> smf(asmeasure(Normal(0.2, σ)), Inf) + csmf(asmeasure(Normal(0.2, σ)), -Inf) + logsmf(asmeasure(Normal(0.2, σ)), Inf)
    @test only(Zygote.gradient(g, 1.3)) == 0 && mooncake_derivative(g, 1.3) == 0
end

@testset "derivatives of truncated gamma and beta distributions" begin
    # Via the differentiable incomplete gamma and beta functions, in the
    # upper tails via their complements; dx/du = mass / pdf(x):
    for (d, a, b) in [(Beta(2.0, 3.0), 0.1, 0.8), (Gamma(2.0, 1.5), 1.0, 4.0), (Gamma(2.0, 1.5), 40.0, Inf)]
        f = transport_to(truncate(asmeasure(d), a .. b), StdUniform())
        @test ForwardDiff.derivative(f, 0.4) ≈ (ccdf(d, a) - ccdf(d, b)) / pdf(d, f(0.4)) rtol = 1e-8
    end
    @test ForwardDiff.derivative(transport_to(truncate(asmeasure(Beta(2.0, 3.0)), 0.1 .. 0.8), StdUniform()), 0.4) ≈ 0.5179695278

    # With respect to the shapes, against central differences:
    mooncake_derivative(f, x) = Mooncake.value_and_gradient!!(Mooncake.prepare_gradient_cache(f, x), f, x)[2][2]
    for (f, x, backends) in [
        (α -> massof(truncate(asmeasure(Gamma(α, 1.5)), 1.0 .. 4.0)), 2.0, (:ForwardDiff,)),
        (α -> massof(truncate(asmeasure(Gamma(α, 1.5)), 40.0 .. Inf)), 2.0, (:ForwardDiff, :Zygote, :Mooncake)),
        (α -> massof(truncate(asmeasure(Beta(α, 3.0)), 0.995 .. 1.0)), 2.0, (:ForwardDiff,)),
        (α -> transport_to(truncate(asmeasure(Gamma(α, 1.5)), 40.0 .. Inf), StdUniform())(0.4), 2.0, (:ForwardDiff,)),
        (α -> transport_to(truncate(asmeasure(Beta(α, 3.0)), 0.1 .. 0.8), StdUniform())(0.4), 2.0, (:ForwardDiff, :Zygote, :Mooncake)),
    ]
        ∂f_ref = (f(x + 1e-6) - f(x - 1e-6)) / 2e-6
        :ForwardDiff in backends && @test ForwardDiff.derivative(f, x) ≈ ∂f_ref rtol = 1e-6
        :Zygote in backends && @test only(Zygote.gradient(f, x)) ≈ ∂f_ref rtol = 1e-6
        :Mooncake in backends && @test mooncake_derivative(f, x) ≈ ∂f_ref rtol = 1e-6
    end
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

# An RNG whose uniform draws are all exactly zero:
struct _AllZeroRNG <: MeasureBase.Random.AbstractRNG end
MeasureBase.Random.rand(::_AllZeroRNG, ::MeasureBase.Random.SamplerTrivial{MeasureBase.Random.CloseOpen01{Float64}}) = 0.0
MeasureBase.Random.rand(::_AllZeroRNG, ::MeasureBase.Random.SamplerTrivial{MeasureBase.Random.CloseOpen01{Float32}}) = 0.0f0

@testset "transported draws of exact zero uniform draws" begin
    using MeasureBase: GenContext, StdLogistic, StdNormal, StdUniform
    # A uniform draw of zero stands for its grid cell (half the smallest
    # nonzero draw), not for floatmin, also in batches:
    for T in (Float32, Float64)
        ctx = GenContext{T}(_AllZeroRNG())
        @test rand(ctx, StdUniform()) == 0
        u = eps(T) / 4
        @test rand(ctx, StdLogistic()) ≈ log(u / (1 - u))
        @test all(≈(log(u / (1 - u))), rand(ctx, StdLogistic()^3))
        t = truncate(StdNormal(), -Inf .. 2.0)
        @test rand(ctx, t) ≈ quantile(Normal(), u * cdf(Normal(), 2.0)) rtol = sqrt(eps(T))
        @test all(isfinite, rand(ctx, t^3))
    end
end
