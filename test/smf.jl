using Distributions: Normal, Logistic, Exponential, Uniform, Poisson, Bernoulli, Cauchy
using Distributions: pdf, cdf, ccdf, logcdf, logccdf, cquantile, invlogcdf, invlogccdf
using IntervalSets
using SpecialFunctions: erfc, erfcinv
using InverseFunctions: inverse

smf_measures = [
    StdNormal()
    StdLogistic()
    StdUniform()
    Half(StdNormal())
    Half(StdLogistic())
    Half(2 * StdNormal())
]

@testset "smf" begin
    for μ in smf_measures
        test_smf(μ)
    end
end

@testset "complementary and log Stieltjes measure functions" begin
    # Against the cdfs of the corresponding distributions, deep into both
    # tails:
    for (μ, d, xs, ps, lps) in (
        (StdNormal(), Normal(), (-40.0, -9.0, -1.0, 0.3, 9.0, 40.0), (1e-300, 1e-20, 0.3, 0.9), (-700.0, -40.0, -1.0, -0.01)),
        (StdLogistic(), Logistic(), (-800.0, -40.0, -1.0, 0.3, 40.0, 800.0), (1e-300, 1e-20, 0.3, 0.9), (-700.0, -40.0, -1.0, -1e-10)),
        (StdExponential(), Exponential(), (0.0, 1e-20, 0.3, 40.0, 800.0), (1e-300, 1e-20, 0.3, 0.9), (-700.0, -40.0, -1.0, -1e-10)),
        (StdUniform(), Uniform(), (0.0, 1e-20, 0.3, 1.0), (0.1, 0.3, 0.9), (-2.0, -1.0, -0.1)),
    )
        for x in xs
            @test smf(μ, x) ≈ cdf(d, x)
            @test csmf(μ, x) ≈ ccdf(d, x)
            @test logsmf(μ, x) ≈ logcdf(d, x)
            @test logcsmf(μ, x) ≈ logccdf(d, x)
        end
        for p in ps
            @test invcsmf(μ, p) ≈ cquantile(d, p)
            @test csmf(μ, invcsmf(μ, p)) ≈ p
        end
        for lp in lps
            @test invlogcsmf(μ, lp) ≈ invlogccdf(d, lp)
            @test invlogsmf(μ, lp) ≈ invlogcdf(d, lp)
        end
        @test isnan(invcsmf(μ, 1.5)) && isnan(invlogcsmf(μ, 0.5)) && isnan(invlogsmf(μ, 0.5))
        @test csmf(μ, 0.3f0) isa Float32 && logcsmf(μ, 0.3f0) isa Float32
    end

    # Weights scale, Half doubles the upper tail of its parent:
    w = 3 * StdNormal()
    @test csmf(w, 9.0) ≈ 3 * ccdf(Normal(), 9.0)
    @test logcsmf(w, 30.0) ≈ log(3) + logccdf(Normal(), 30.0)
    @test invcsmf(w, 3e-20) ≈ cquantile(Normal(), 1e-20)
    @test invlogcsmf(w, log(3) - 40) ≈ invlogccdf(Normal(), -40.0)
    @test csmf(Half(StdNormal()), 9.0) ≈ 2 * ccdf(Normal(), 9.0)
    @test invcsmf(Half(StdNormal()), 2e-20) ≈ cquantile(Normal(), 1e-20)
    # and has the mass of its parent:
    h = Half(2 * StdNormal())
    @test smf(h, 0.3) ≈ 2 * (2 * cdf(Normal(), 0.3) - 1) && smf(h, Inf) ≈ 2
    @test csmf(h, 9.0) ≈ 4 * ccdf(Normal(), 9.0)
    @test invsmf(h, 1.0) ≈ cquantile(Normal(), 0.25) && invcsmf(h, 4e-20) ≈ cquantile(Normal(), 1e-20)
    @test isnan(invsmf(h, 2.5)) && isnan(invcsmf(h, -0.1))
    @test csmf(Lebesgue(0.0 .. 2.0), 0.5) == 1.5 && invcsmf(Lebesgue(0.0 .. 2.0), 1.5) == 0.5
    @test csmf(Lebesgue(), 0.5) isa MeasureBase.NoSMF

    # Masses of intervals in the upper tail:
    @test massof(StdNormal(), 9.0 .. Inf) ≈ ccdf(Normal(), 9.0)
    @test massof(StdNormal(), 9.0 .. 10.0) ≈ ccdf(Normal(), 9.0) - ccdf(Normal(), 10.0)
    @test massof(StdNormal(), -10.0 .. -9.0) ≈ cdf(Normal(), -9.0) - cdf(Normal(), -10.0)
end

@testset "log-space inverses of the StdNormal Stieltjes measure functions" begin
    # Beyond the underflow of the cdf, Float64 and Float32:
    for lp in (-1e-20, -0.1, -1.0, -40.0, -700.0, -745.0, -800.0, -5000.0, -1e6)
        x = invlogsmf(StdNormal(), lp)
        @test logsmf(StdNormal(), x) ≈ lp
        @test invlogcsmf(StdNormal(), lp) == -x
    end
    for lp in (-1.0, -40.0, -700.0, -800.0)
        @test invlogsmf(StdNormal(), lp) ≈ invlogcdf(Normal(), lp)
    end
    @test invlogsmf(StdNormal(), logsmf(StdNormal(), -40.0)) ≈ -40
    for lp in (-1.0f0, -40.0f0, -87.0f0, -90.0f0, -200.0f0)
        x = invlogsmf(StdNormal(), lp)
        @test x isa Float32
        @test x ≈ invlogcdf(Normal(), Float64(lp)) rtol = 1e-6
    end
    @test isnan(invlogsmf(StdNormal(), 0.1)) && invlogsmf(StdNormal(), -Inf) == -Inf
    # Inherited by weighted measures and Half:
    @test invlogsmf(2 * StdNormal(), log(2) - 800) ≈ invlogcdf(Normal(), -800.0)
    @test logcsmf(Half(StdNormal()), 40.0) ≈ log(2) + logccdf(Normal(), 40.0)
    @test invlogcsmf(Half(StdNormal()), log(2) - 800) ≈ invlogccdf(Normal(), -800.0)
    @test isnan(invlogcsmf(Half(StdNormal()), 0.1))
    # Conversions between standard measures in the upper tails:
    @test MeasureBase.stdconvert(StdExponential, StdNormal, 13.0f0) ≈ -logccdf(Normal(), 13.0)
    @test MeasureBase.stdconvert(StdNormal, StdExponential, 800.0) ≈ invlogccdf(Normal(), -800.0)
    @test isfinite(MeasureBase.stdconvert(StdExponential, StdNormal, 1e200))
end

@testset "interval masses of measures with atoms" begin
    # Stieltjes measure functions give the masses of left-open, right-closed
    # intervals:
    p = asmeasure(Poisson(2.0))
    @test massof(p, Interval{:open,:closed}(1, 5)) ≈ sum(k -> pdf(Poisson(2.0), k), 2:5)
    @test massof(p, Interval{:open,:closed}(-1, 0)) ≈ pdf(Poisson(2.0), 0)
    @test massof(asmeasure(Bernoulli(0.3)), Interval{:open,:closed}(-1, 0)) ≈ 0.7
    # Other intervals would lose or gain the atoms at their ends:
    @test_throws ArgumentError massof(p, 2 .. 5)
    @test_throws ArgumentError p(0 .. 0)
    @test_throws ArgumentError massof(asmeasure(Bernoulli(0.3)), 0 .. 1)
    @test massof(StdNormal(), 0 .. 1) ≈ cdf(Normal(), 1) - cdf(Normal(), 0)
end

@testset "interval masses" begin
    # Against references in extended precision, in the bulk and the tails,
    # for intervals far shorter than the resolution of the Stieltjes measure
    # functions and in both directions. The tolerances allow for the
    # conditioning of the densities in the tails (and the references need
    # the precision to resolve the shortest interval):
    big_smfs = (
        (StdNormal(), x -> erfc(-big(x) / sqrt(big(2))) / 2, x -> exp(-x^2 / 2) / sqrt(2π)),
        (StdLogistic(), x -> 1 / (1 + exp(-big(x))), x -> exp(-abs(x)) / (1 + exp(-abs(x)))^2),
        (StdExponential(), x -> -expm1(-max(big(x), 0)), x -> exp(-x)),
    )
    setprecision(BigFloat, 1200) do
        for (μ, F, f) in big_smfs
            for a in (-8.0, -1.0, -0.3, 0.0, 0.5, 3.0, 8.0), h in (1e-300, 1e-17, 1e-12, 1e-6, 0.01, 0.1, 1.0, 5.0, -1e-12, -0.3)
                x = a + h
                (μ isa StdExponential && min(a, x) < 0) && continue
                p = smf(μ, a, x)
                @test p ≈ F(x) - F(a) rtol = 16 * eps() * max(1, x^2)
                # Inverses within the resolution of `x` given `p`:
                Δx = abs(p) / f(x)
                Δx < 1 && @test abs(invsmf(μ, a, p) - x) <= 4 * eps() * (abs(x) + Δx)
            end
        end
    end

    for μ in (StdNormal(), StdLogistic(), StdExponential())
        @test smf(μ, 0.5, 0.5) == 0 && smf(μ, Inf, Inf) == 0
        @test smf(μ, 1.5, 0.5) ≈ -smf(μ, 0.5, 1.5)
        @test smf(μ, 0.5, Inf) ≈ csmf(μ, 0.5) && smf(μ, -Inf, 0.5) ≈ smf(μ, 0.5)
        @test smf(μ, -Inf, Inf) == 1
        # NaN outside of the range of the interval masses, finite at its
        # ends:
        @test all(isnan, invsmf.(μ, 0.5, [-1.1 * smf(μ, 0.5), 1.1 * csmf(μ, 0.5), NaN]))
        @test isfinite(invsmf(μ, 0.5, csmf(μ, 0.5))) && invsmf(μ, 0.5, csmf(μ, 0.5)) > invsmf(μ, 0.5, 0.999 * csmf(μ, 0.5))
        @test isfinite(invsmf(μ, -Inf, 0.0)) && isfinite(invsmf(μ, Inf, 0.0))
        @test invsmf(μ, -Inf, smf(μ, 0.5)) ≈ 0.5 && invsmf(μ, Inf, -csmf(μ, 0.5)) ≈ 0.5
        @test smf(μ, 0.5f0, 0.75f0) isa Float32 && invsmf(μ, 0.5f0, 0.1f0) isa Float32
        @test smf(μ, 0, 1) ≈ smf(μ, 0.0, 1.0)
    end
    @test -40 < invsmf(StdNormal(), 0.5, -smf(StdNormal(), 0.5)) < -37

    # Weights and normalization scale them, other measures use the
    # difference of their Stieltjes measure functions:
    @test smf(2 * StdNormal(), 0.0, 1e-12) ≈ 2 * smf(StdNormal(), 0.0, 1e-12)
    @test invsmf(2 * StdNormal(), 0.0, 2e-12) ≈ invsmf(StdNormal(), 0.0, 1e-12)
    n = normalize(2 * Half(StdNormal()))
    @test smf(n, 0.0, 1e-12) ≈ smf(Half(StdNormal()), 0.0, 1e-12) ≈ 2 * smf(StdNormal(), 0.0, 1e-12)
    @test invsmf(n, 0.0, smf(n, 0.0, 1e-12)) ≈ 1e-12
    c = asmeasure(Cauchy())
    @test smf(c, 3.0, 4.0) ≈ cdf(Cauchy(), 4.0) - cdf(Cauchy(), 3.0)
    @test smf(c, 30.0, 40.0) ≈ ccdf(Cauchy(), 30.0) - ccdf(Cauchy(), 40.0)
    @test invsmf(c, 3.0, smf(c, 3.0, 4.0)) ≈ 4.0 && invsmf(c, 30.0, smf(c, 30.0, 40.0)) ≈ 40.0
    @test smf(Lebesgue(), 0.5, 2.0) == 1.5 && invsmf(Lebesgue(), 0.5, 1.5) == 2.0
    @test smf(Lebesgue(0.0 .. 2.0), 0.5, 3.0) == 1.5
    @test smf(StdNormal()^2, 0.0, 1.0) isa MeasureBase.NoSMF
    @test invsmf(StdNormal()^2, 0.0, 0.5) isa MeasureBase.NoSMFInverse
end

@testset "Half and truncations close to their lower end" begin
    # Transports of variates just above the lower end, against references
    # in extended precision:
    Φ(x) = erfc(-big(x) / sqrt(big(2))) / 2
    Φinv(p) = -sqrt(big(2)) * erfcinv(2 * big(p))
    setprecision(BigFloat, 1200) do
        for (μ, a, b) in ((Half(StdNormal()), 0.0, Inf), (truncate(StdNormal(), 0.0 .. Inf), 0.0, Inf), (truncate(StdNormal(), 0.5 .. 2.0), 0.5, 2.0))
            f = transport_to(StdNormal(), μ)
            for h in (1e-300, 1e-17, 1e-12, 1e-6)
                x = a + h
                x == a && continue
                z = f(x)
                @test z ≈ Φinv((Φ(x) - Φ(a)) / (Φ(b) - Φ(a))) rtol = 1e-14
                # The round trip, limited by the conditioning of the uniform pivot:
                @test inverse(f)(z) ≈ x rtol = 1e-12
            end
        end
    end
    h = Half(StdNormal())
    for x in (1e-300, 1e-17, 1e-10, 0.3)
        @test invsmf(h, smf(h, x)) ≈ x rtol = 4 * eps()
    end
end
