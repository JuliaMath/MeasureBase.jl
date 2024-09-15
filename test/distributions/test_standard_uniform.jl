# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using DistributionMeasures
using Test

using Random, Statistics, LinearAlgebra
using Distributions, PDMats
using StableRNGs
using FillArrays
using ForwardDiff


@testset "StandardDist{Uniform}" begin
    stblrng() = StableRNG(789990641)

    @testset "StandardDist{Uniform,0}" begin
        @test @inferred(Uniform(StandardDist{Uniform}())) isa Uniform{Float64}
        @test @inferred(Uniform(StandardDist{Uniform}())) == Uniform()
        @test @inferred(convert(Uniform, StandardDist{Uniform}())) == Uniform()

        d = StandardDist{Uniform}()
        dref = Uniform()

        @test @inferred(minimum(d)) == minimum(dref)
        @test @inferred(maximum(d)) == maximum(dref)
        
        @test @inferred(Distributions.params(d)) == ()
        @test @inferred(partype(d)) == partype(dref)
        
        @test @inferred(location(d)) == location(dref)
        @test @inferred(scale(d)) == scale(dref)
        
        @test @inferred(eltype(typeof(d))) == eltype(typeof(dref))
        @test @inferred(eltype(d)) == eltype(dref)

        @test @inferred(length(d)) == length(dref)
        @test @inferred(size(d)) == size(dref)

        @test @inferred(mean(d)) == mean(dref)
        @test @inferred(median(d)) == median(dref)
        @test @inferred(mode(d)) == mode(dref)
        @test @inferred(modes(d)) ≈ modes(dref)
        
        @test @inferred(var(d)) ≈ var(dref)
        @test @inferred(std(d)) ≈ std(dref)
        @test @inferred(skewness(d)) == skewness(dref)
        @test @inferred(kurtosis(d)) ≈ kurtosis(dref)
        
        @test @inferred(entropy(d)) == entropy(dref)
        
        for x in [-0.5, 0.0, 0.25, 0.75, 1.0, 1.5]
            @test @inferred(logpdf(d, x)) == logpdf(dref, x)
            @test @inferred(pdf(d, x)) == pdf(dref, x)
            @test @inferred(logcdf(d, x)) == logcdf(dref, x)
            @test @inferred(cdf(d, x)) == cdf(dref, x)
            @test @inferred(logccdf(d, x)) == logccdf(dref, x)
            @test @inferred(ccdf(d, x)) == ccdf(dref, x)
        end

        for p in [0.0, 0.25, 0.75, 1.0]
            @test @inferred(quantile(d, p)) == quantile(dref, p)
            @test @inferred(cquantile(d, p)) == cquantile(dref, p)
        end

        for t in [-3, 0, 3]
            @test @inferred(mgf(d, t)) == mgf(dref, t)
            @test @inferred(cf(d, t)) == cf(dref, t)
        end

        @test @inferred(rand(stblrng(), d)) == rand(stblrng(), dref)
        @test @inferred(rand!(stblrng(), d, fill(0.0))) == rand!(stblrng(), dref, fill(0.0))
        @test @inferred(rand(stblrng(), d, 5)) == rand(stblrng(), dref, 5)

        @test @inferred(truncated(StandardDist{Uniform}(), -0.5f0, 0.7f0)) isa Uniform{Float64}
        @test truncated(StandardDist{Uniform}(), -0.5f0, 0.7f0) == Uniform(0.0f0, 0.7f0)
        @test truncated(StandardDist{Uniform}(), 0.2f0, 0.7f0) == Uniform(0.2f0, 0.7f0)

        @test @inferred(product_distribution(fill(StandardDist{Uniform}(), 3))) isa DistributionMeasures.StandardDist{Uniform,1}
        @test product_distribution(fill(StandardDist{Uniform}(), 3)) == DistributionMeasures.StandardDist{Uniform}(3)
    end


    @testset "StandardDist{Uniform,1}" begin
        d = DistributionMeasures.StandardDist{Uniform}(3)
        dref = product_distribution(fill(Uniform(), 3))

        @test @inferred(eltype(typeof(d))) == eltype(typeof(dref))
        @test @inferred(eltype(d)) == eltype(dref)

        @test @inferred(length(d)) == length(dref)
        @test @inferred(size(d)) == size(dref)

        @test @inferred(Distributions.params(d)) == ()
        @test @inferred(partype(d)) == partype(dref)

        @test @inferred(mean(d)) == mean(dref)
        @test @inferred(var(d)) ≈ var(dref)
        @test @inferred(cov(d)) ≈ cov(dref)

        @test @inferred(mode(d)) == [0.5, 0.5, 0.5]
        @test @inferred(modes(d)) == fill([0, 0,0 ])

        @test @inferred(invcov(d)) == inv(cov(dref))
        @test @inferred(logdetcov(d)) == logdet(cov(dref))

        @test @inferred(entropy(d)) == entropy(dref)

        for x in fill.([-Inf, -1.3, 0.0, 1.3, +Inf], 3)
            @test @inferred(Distributions.insupport(d, x)) == Distributions.insupport(dref, x)
            @test @inferred(logpdf(d, x)) == logpdf(dref, x)
            @test @inferred(pdf(d, x)) == pdf(dref, x)
            @test @inferred(gradlogpdf(d, x)) == ForwardDiff.gradient(x -> logpdf(d, x), x)
        end

        @test @inferred(rand(stblrng(), d)) == rand(stblrng(), d)
        @test @inferred(rand!(stblrng(), d, zeros(3))) == rand!(stblrng(), d, zeros(3))
        @test @inferred(rand!(stblrng(), d, zeros(3, 10))) == rand!(stblrng(), d, zeros(3, 10))
    end
end
