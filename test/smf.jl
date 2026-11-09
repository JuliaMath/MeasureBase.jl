using Distributions: Normal, Logistic, Exponential, Uniform, Poisson, Bernoulli
using Distributions: pdf, cdf, ccdf, logcdf, logccdf, cquantile, invlogcdf, invlogccdf
using IntervalSets

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
