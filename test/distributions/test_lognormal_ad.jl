# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test, MeasureBase, Distributions
using MeasureBase: normalize
using IntervalSets
import ForwardDiff, Enzyme

@testset "LogNormal support derivatives" begin
    parameters = [log(4.0), log(0.7)]
    for f in (smf, csmf, logsmf, logcsmf), endpoint in (0.0, Inf)
        g(p) = f(asmeasure(LogNormal(p[1], exp(p[2]))), endpoint)
        @test ForwardDiff.gradient(g, parameters) == [0.0, 0.0]
    end
    law = asmeasure(LogNormal(parameters[1], exp(parameters[2])))
    @test ForwardDiff.derivative(x -> smf(law, x), 0.0) == 0.0
    for (f, expected) in ((smf, 1.0), (csmf, 0.0), (logsmf, 0.0), (logcsmf, -Inf))
        @test f(asmeasure(LogNormal(0.0, 0.0)), 1.0) == expected
    end
    for x in (0, 1)
        @test (@inferred smf(asmeasure(LogNormal(0.0f0, 1.0f0)), x)) ≈ cdf(LogNormal(), x)
    end

    probability = 0.5
    sigma = exp(parameters[2])
    z = (log(8.0)-parameters[1])/sigma
    q = quantile(Normal(), probability*cdf(Normal(), z))
    value = exp(parameters[1]+sigma*q)
    ratio = probability*pdf(Normal(), z)/pdf(Normal(), q)
    expected = value .* [1-ratio, sigma*(q-z*ratio)]
    f(p) = invsmf(
        normalize(truncate(asmeasure(LogNormal(p[1], exp(p[2]))), 0.0 .. 8.0)),
        probability,
    )
    @test f(parameters) ≈ value
    @test ForwardDiff.gradient(f, parameters) ≈ expected rtol=1e-12
    @test only(Enzyme.gradient(Enzyme.Reverse, f, parameters)) ≈ expected rtol=1e-12
end
