# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

import Zygote
import ForwardDiff

using MeasureBase
using MeasureBase: transport_to, StdUniform, StdNormal, StdExponential, StdLogistic
using MeasureBase: mintegrate_exp, logdensities


# Log-densities of the standard measures, as independent references:
_ref_logpdf(::StdNormal, x) = -x^2 / 2 - log(2π) / 2
_ref_logpdf(::StdExponential, x) = -x
_ref_logpdf(::StdLogistic, x) = -abs(x) - 2 * log1p(exp(-abs(x)))

@testset "Zygote" begin
    @testset "standard conversions" begin
        # The derivative of a transport is the ratio of the source density at
        # the variate and the target density at its image. Includes the
        # branch points and the tails, where the branches not taken are
        # singular:
        for (ν, μ, xs) in [
            (StdExponential(), StdNormal(), [-40.0, -8.5, -1.0, -0.0, 0.0, 1.0, 8.5, 9.0, 20.0, 37.0]),
            (StdNormal(), StdExponential(), [1e-20, 1e-5, 0.5, prevfloat(log(2.0)), log(2.0), 1.0, 9.0, 40.0, 700.0, 800.0]),
            (StdNormal(), StdLogistic(), [-700.0, -40.0, -1.0, -0.0, 0.0, 1.0, 40.0, 700.0]),
        ]
            f = transport_to(ν, μ)
            for x in xs
                ref = exp(_ref_logpdf(μ, x) - _ref_logpdf(ν, f(x)))
                @test Zygote.gradient(f, x)[1] ≈ ref rtol = 1e-10
                @test ForwardDiff.derivative(f, x) ≈ ref rtol = 1e-10
            end
        end
        # Saturated beyond floatmin tail probabilities:
        @test Zygote.gradient(transport_to(StdNormal(), StdLogistic()), 800.0)[1] == 0
        L = [-9.0, 0.0, 9.0, 40.0]
        f = transport_to(StdNormal()^4, StdLogistic()^4)
        @test Zygote.gradient(L -> sum(f(L)), L)[1] ≈ exp.(_ref_logpdf.(StdLogistic(), L) .- _ref_logpdf.(StdNormal(), f(L)))
    end

    # The batched kernels must not accumulate into result buffers, which
    # Zygote can't differentiate and which can't hold dual numbers of
    # captured parameters:
    @testset "batched densities with out-of-support variates" begin
        # log-density -|x|² - sum(x) inside the support, -Inf outside:
        post = mintegrate_exp(x -> -sum(abs2, x), StdExponential()^2)
        X = [0.1 -0.5 1.0; 0.3 0.2 2.0]
        ref = [-1.2 0.0 -3.0; -1.6 0.0 -5.0]
        f = X -> sum(logdensities(post, X)[[1, 3]])
        @test ForwardDiff.gradient(f, X) ≈ ref
        @test Zygote.gradient(f, X)[1] ≈ ref

        # Of the variates 0.5 and 2.0 only the first is inside the support:
        g = a -> sum(logdensities(mintegrate_exp(x -> -a * x^2, StdUniform()), [0.5, 2.0]))
        @test ForwardDiff.derivative(g, 1.0) ≈ -0.25
        @test Zygote.gradient(g, 1.0)[1] ≈ -0.25
    end
end


# Wrapped distributions under Zygote without ForwardDiffPullbacks (which the
# Distributions tests load later), with the distributions built inside the
# differentiated functions, against central differences of Distributions'
# own functions:
import Distributions
using Distributions: Chisq, FDist, InverseGamma, BetaPrime, Normal, MixtureModel
using IntervalSets
using MeasureBase: asmeasure, normalize

@testset "Zygote with Distributions" begin
    @test isnothing(Base.get_extension(MeasureBase, :MeasureBaseDistributionsForwardDiffPullbacksExt))
    fd(f, a) = (h = 1e-4 * a; (8 * (f(a + h) - f(a - h)) - (f(a + 2h) - f(a - 2h))) / (12 * h))
    to_normal(d, x) = Distributions.quantile(Normal(), Distributions.cdf(d, x))
    from_normal(d, z) = Distributions.quantile(d, Distributions.cdf(Normal(), z))
    mix(w) = MixtureModel([Normal(-1.0, 1.0), Normal(2.0, 3.0)], [w, 1 - w])
    for (f, f_ref, a) in [
        (ν -> transport_to(StdNormal(), asmeasure(Chisq(ν)))(0.7), ν -> to_normal(Chisq(ν), 0.7), 2.3),
        (x -> transport_to(StdNormal(), asmeasure(Chisq(2.3)))(x), x -> to_normal(Chisq(2.3), x), 0.7),
        (ν -> transport_to(StdNormal(), asmeasure(FDist(ν, 5.0)))(0.8), ν -> to_normal(FDist(ν, 5.0), 0.8), 3.0),
        (ν -> transport_to(asmeasure(FDist(ν, 5.0)), StdNormal())(0.3), ν -> from_normal(FDist(ν, 5.0), 0.3), 3.0),
        (α -> transport_to(StdNormal(), asmeasure(InverseGamma(α, 1.3)))(0.8), α -> to_normal(InverseGamma(α, 1.3), 0.8), 2.5),
        (β -> transport_to(asmeasure(BetaPrime(2.0, β)), StdNormal())(0.3), β -> from_normal(BetaPrime(2.0, β), 0.3), 3.0),
        (ν -> logdensityof(normalize(truncate(asmeasure(Chisq(ν)), 0.5 .. 3.0)), 1.0),
            ν -> Distributions.logpdf(Distributions.truncated(Chisq(ν), 0.5, 3.0), 1.0), 2.3),
        # Mixture quantiles are roots, their bisection isn't differentiated:
        (z -> transport_to(asmeasure(mix(0.3)), StdNormal())(z), z -> from_normal(mix(0.3), z), 0.4),
        (w -> transport_to(asmeasure(mix(w)), StdNormal())(0.4), w -> from_normal(mix(w), 0.4), 0.3),
    ]
        @test Zygote.gradient(f, a)[1] ≈ fd(f_ref, a) rtol = 1e-8
    end
end
