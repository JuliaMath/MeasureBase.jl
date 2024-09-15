# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using DistributionMeasures
using Test

using Random, Statistics, LinearAlgebra
using Distributions, PDMats
using StableRNGs


@testset "StandardDist{Normal}" begin
    stblrng() = StableRNG(789990641)

    @testset "StandardDist{Normal,0}" begin
        @test @inferred(Normal(StandardDist{Normal}())) isa Normal{Float64}
        @test @inferred(Normal(StandardDist{Normal}())) == Normal()
        @test @inferred(convert(Normal, StandardDist{Normal}())) == Normal()

        d = StandardDist{Normal}()
        dref = Normal()

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
        
        @test @inferred(var(d)) == var(dref)
        @test @inferred(std(d)) == std(dref)
        @test @inferred(skewness(d)) == skewness(dref)
        @test @inferred(kurtosis(d)) == kurtosis(dref)
        
        @test @inferred(entropy(d)) == entropy(dref)
        
        for x in [-Inf, -1.3, 0.0, 1.3, +Inf]
            @test @inferred(gradlogpdf(d, x)) == gradlogpdf(dref, x)
            
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

        @test @inferred(truncated(StandardDist{Normal}(), -2.2f0, 3.1f0)) isa Truncated{Normal{Float64}}
        @test truncated(StandardDist{Normal}(), -2.2f0, 3.1f0) == truncated(Normal(0.0, 1.0), -2.2f0, 3.1f0)

        @test @inferred(product_distribution(fill(StandardDist{Normal}(), 3))) isa StandardDist{Normal,1}
        @test product_distribution(fill(StandardDist{Normal}(), 3)) == StandardDist{Normal}(3)
    end


    @testset "StandardDist{Normal,1}" begin
        @test @inferred(StandardDist{Normal}(3)) isa StandardDist{Normal,1}
        @test @inferred(StandardDist{Normal}(3)) isa StandardDist{Normal,1}
        @test @inferred(StandardDist{Normal}(3)) isa StandardDist{Normal,1}

        @test @inferred(MvNormal(StandardDist{Normal}(3))) isa MvNormal{Int}
        @test @inferred(MvNormal(StandardDist{Normal}(3))) == MvNormal(ScalMat(3, 1.0))
        @test @inferred(convert(MvNormal, StandardDist{Normal}(3))) == MvNormal(ScalMat(3, 1.0))

        d = StandardDist{Normal}(3)
        dref = MvNormal(ScalMat(3, 1.0))

        @test @inferred(eltype(typeof(d))) == eltype(typeof(dref))
        @test @inferred(eltype(d)) == eltype(dref)

        @test @inferred(length(d)) == length(dref)
        @test @inferred(size(d)) == size(dref)

        @test @inferred(Distributions.params(d)) == ()
        @test @inferred(partype(d)) == partype(dref)

        @test @inferred(mean(d)) == mean(dref)
        @test @inferred(var(d)) == var(dref)
        @test @inferred(cov(d)) == cov(dref)
        
        @test @inferred(mode(d)) == mode(dref)
        @test @inferred(modes(d)) == modes(dref)

        @test @inferred(invcov(d)) == invcov(dref)
        @test @inferred(logdetcov(d)) == logdetcov(dref)

        @test @inferred(entropy(d)) == entropy(dref)

        for x in fill.([-Inf, -1.3, 0.0, 1.3, +Inf], 3)
            # Distributions.insupport is inconsistent at +- Inf between Normal and MvNormal
            if !any(isinf, x)
                @test @inferred(Distributions.insupport(d, x)) == Distributions.insupport(dref, x)
            end
            @test @inferred(logpdf(d, x)) == logpdf(dref, x)
            @test @inferred(pdf(d, x)) == pdf(dref, x)
            @test @inferred(sqmahal(d, x)) == sqmahal(dref, x)
            @test @inferred(gradlogpdf(d, x)) == gradlogpdf(dref, x)
        end

        @test @inferred(rand(stblrng(), d)) == rand(stblrng(), d)
        @test @inferred(rand!(stblrng(), d, zeros(3))) == rand!(stblrng(), d, zeros(3))
        @test @inferred(rand!(stblrng(), d, zeros(3, 10))) == rand!(stblrng(), d, zeros(3, 10))
    end
end
