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

    # The batched kernels accumulate into result buffers, which Zygote can't
    # differentiate and which can't hold dual numbers of captured parameters
    # (marked broken until the kernels allocate their results):
    @testset "batched densities with out-of-support variates" begin
        # log-density -|x|² - sum(x) inside the support, -Inf outside:
        post = mintegrate_exp(x -> -sum(abs2, x), StdExponential()^2)
        X = [0.1 -0.5 1.0; 0.3 0.2 2.0]
        ref = [-1.2 0.0 -3.0; -1.6 0.0 -5.0]
        f = X -> sum(logdensities(post, X)[[1, 3]])
        @test ForwardDiff.gradient(f, X) ≈ ref
        @test Zygote.gradient(f, X)[1] ≈ ref broken = true

        # Of the variates 0.5 and 2.0 only the first is inside the support:
        g = a -> sum(logdensities(mintegrate_exp(x -> -a * x^2, StdUniform()), [0.5, 2.0]))
        @test ForwardDiff.derivative(g, 1.0) ≈ -0.25 broken = true
        @test Zygote.gradient(g, 1.0)[1] ≈ -0.25 broken = true
    end
end
