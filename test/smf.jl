using Distributions: Normal, Logistic, Exponential, Uniform
using Distributions: cdf, ccdf, logcdf, logccdf, cquantile, invlogcdf, invlogccdf
using IntervalSets

smf_measures = [
    StdNormal()
    StdLogistic()
    StdUniform()
    Half(StdNormal())
    Half(StdLogistic())
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
    @test csmf(Lebesgue(0.0 .. 2.0), 0.5) == 1.5 && invcsmf(Lebesgue(0.0 .. 2.0), 1.5) == 0.5
    @test csmf(Lebesgue(), 0.5) isa MeasureBase.NoSMF

    # Masses of intervals in the upper tail:
    @test massof(StdNormal(), 9.0 .. Inf) ≈ ccdf(Normal(), 9.0)
    @test massof(StdNormal(), 9.0 .. 10.0) ≈ ccdf(Normal(), 9.0) - ccdf(Normal(), 10.0)
    @test massof(StdNormal(), -10.0 .. -9.0) ≈ cdf(Normal(), -9.0) - cdf(Normal(), -10.0)
end
