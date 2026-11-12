# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: normalize
using Distributions
using LinearAlgebra: Diagonal
using AffineMaps: MulAdd, Mul
using IntervalSets
import ForwardDiff, Zygote, Mooncake

@testset "moments of distributions" begin
    for d in [Normal(0.3, 1.7), Gamma(2.0, 1.5), Beta(2.0, 3.0), Poisson(3.5), Laplace(0.2, 0.4)]
        m = asmeasure(d)
        @test mean(m) ≈ mean(d)
        @test var(m) ≈ var(d)
        @test std(m) ≈ std(d)
        @test mean(2.5 * m) ≈ mean(d) && var(normalize(2.5 * m)) ≈ var(d)
        @test_throws ArgumentError cov(m)
    end

    for d in [MvNormal([0.3, -0.2], [2.0 0.5; 0.5 1.0]), Dirichlet([1.0, 2.0, 3.0])]
        m = asmeasure(d)
        @test mean(m) ≈ mean(d)
        @test var(m) ≈ var(d)
        @test cov(m) ≈ cov(d)
    end

    ds = [Normal(0.3, 1.7), Gamma(2.0, 1.5), Exponential(0.5)]
    pd = product_distribution(ds)
    m = asmeasure(pd)
    @test m isa MeasureBase.ProductMeasure
    @test mean(m) ≈ mean(pd)
    @test var(m) ≈ var(pd)
    @test cov(m) ≈ cov(pd)
    @test mean(productmeasure((a = asmeasure(ds[1]), b = asmeasure(ds[2])))) == (a = mean(ds[1]), b = mean(ds[2]))
    @test var(asmeasure(Normal(0.3, 1.7))^3) ≈ fill(1.7^2, 3)

    # Affine pushforwards:
    A, b = [2.0 0.5; -1.0 1.5], [0.3, -0.7]
    mvn = MvNormal([0.3, -0.2], [2.0 0.5; 0.5 1.0])
    ν = pushfwd(MulAdd(A, b), asmeasure(mvn))
    @test mean(ν) ≈ A * mean(mvn) + b
    @test cov(ν) ≈ A * cov(mvn) * A'
    @test var(ν) ≈ var(MvNormal(A * mean(mvn) + b, A * cov(mvn) * A'))
    @test mean(pushfwd(Mul(2.0), asmeasure(Gamma(2.0, 1.5)))) ≈ 2 * mean(Gamma(2.0, 1.5))
    @test var(pushfwd(Mul(2.0), asmeasure(Gamma(2.0, 1.5)))) ≈ 4 * var(Gamma(2.0, 1.5))
end

@testset "moments of native truncated normal distributions" begin
    t = truncate(asmeasure(Normal(0, 2)), -3 .. 3)
    v = 4 * (1 - 3 * pdf(Normal(), 1.5) / (2cdf(Normal(), 1.5) - 1))
    @test mean(t) == mean(normalize(t)) == 0
    @test var(t) ≈ v && std(normalize(t)) ≈ sqrt(v)
    t32 = truncate(asmeasure(Normal(0.0f0, 2.0f0)), -3.0f0 .. 3.0f0)
    @test mean(t32) === 0.0f0 && var(t32) ≈ Float32(v)
    narrow32 = truncate(asmeasure(Normal(0.2f0, 1.3f0)), 1.0f0 .. 1.000001f0)
    @test var(narrow32) ≈ (narrow32.hi - narrow32.lo)^2 / 12 rtol = sqrt(eps(Float32))

    # Both parameters, including the constant infinite endpoints and the
    # centered series for a short interval, against central differences.
    for (a, b) in [(-1.0, 2.0), (0.0, Inf), (-Inf, 0.0), (-Inf, Inf), (1.0, 1.02)]
        for moment in (mean, var)
            f(p) = moment(normalize(truncate(asmeasure(Normal(p[1], p[2])), a .. b)))
            p = [0.2, 1.3]
            reference = [(f(p + h) - f(p - h)) / 2e-5 for h in ([1e-5, 0.0], [0.0, 1e-5])]
            @test ForwardDiff.gradient(f, p) ≈ reference rtol = 1e-6 atol = 1e-9
            @test only(Zygote.gradient(f, p)) ≈ reference rtol = 1e-6 atol = 1e-9
            cache = Mooncake.prepare_gradient_cache(f, p)
            @test Mooncake.value_and_gradient!!(cache, f, p)[2][2] ≈ reference rtol = 1e-6 atol =
                1e-9
        end
    end
end
