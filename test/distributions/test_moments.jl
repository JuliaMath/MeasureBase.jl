# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: normalize
using Distributions
using LinearAlgebra: Diagonal
using AffineMaps: MulAdd, Mul

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
