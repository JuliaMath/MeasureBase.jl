# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using LinearAlgebra
using Distributions, ArraysOfArrays
using FillArrays: Fill
import ForwardDiff, Zygote
using MeasureBase: StdNormal, StdUniform, StdLogistic, transport_to


@testset "trafo_utils" begin
    xs = rand(5)
    @test Zygote.jacobian(MeasureBase._pushfront, xs, 42)[1] ≈ ForwardDiff.jacobian(xs -> MeasureBase._pushfront(xs, 1), xs)
    @test Zygote.jacobian(MeasureBase._pushfront, xs, 42)[2] ≈ vec(ForwardDiff.jacobian(x -> MeasureBase._pushfront(xs, x[1]), [42]))
    @test Zygote.jacobian(MeasureBase._pushback, xs, 42)[1] ≈ ForwardDiff.jacobian(xs -> MeasureBase._pushback(xs, 1), xs)
    @test Zygote.jacobian(MeasureBase._pushback, xs, 42)[2] ≈ vec(ForwardDiff.jacobian(x -> MeasureBase._pushback(xs, x[1]), [42]))
    @test Zygote.jacobian(MeasureBase._rev_cumsum, xs)[1] ≈ ForwardDiff.jacobian(MeasureBase._rev_cumsum, xs)
    @test Zygote.jacobian(MeasureBase._exp_cumsum_log, xs)[1] ≈ ForwardDiff.jacobian(MeasureBase._exp_cumsum_log, xs) ≈ ForwardDiff.jacobian(cumprod, xs)
end


@testset "Zygote gradients of densities" begin
    m = asmeasure(Normal())
    @test Zygote.gradient(x -> -logdensityof(m, x), 1.0)[1] ≈ 1.0
    @test Zygote.gradient(x -> 3 * logdensityof(m, x), 1.0)[1] ≈ -3.0
    @test Zygote.gradient(x -> -logdensityof(m^2, x), [1.0, 2.0])[1] ≈ [1.0, 2.0]
    @test Zygote.gradient(x -> -logdensityof(StdNormal(), x), 1.0)[1] ≈ 1.0
    @test Zygote.gradient(x -> -logdensityof(Lebesgue(), x), 1.0)[1] ∈ (0, nothing)
end


@testset "Zygote gradients of univariate transports" begin
    for (ν, μ, x) in [
        (StdNormal(), Normal(1, 2), 0.7), (Normal(1, 2), StdNormal(), 0.7),
        (StdNormal(), Gamma(2.0), 1.3), (Gamma(2.0), StdNormal(), 0.4),
        (StdLogistic(), TDist(3), 0.8), (TDist(3), StdUniform(), 0.3),
        (Beta(2, 3), Weibull(0.7, 1.3), 0.9), (asmeasure(Logistic()), asmeasure(Normal(0.5)), 0.2),
    ]
        f = transport_to(ν, μ)
        @test Zygote.gradient(f, x)[1] ≈ ForwardDiff.derivative(f, x)
    end
end


@testset "Zygote gradients of transports in the tails" begin
    # Through the standard conversions, against the density ratios:
    φ(z) = exp(-z^2 / 2) / sqrt(2π)
    for (ν, μ, x, ref) in [
        (StdNormal(), Exponential(1.0), 40.0, z -> pdf(Exponential(1.0), 40.0) / φ(z)),
        (StdNormal(), Weibull(1.5, 2.0), 30.0, z -> pdf(Weibull(1.5, 2.0), 30.0) / φ(z)),
        (StdNormal(), Gumbel(), 40.0, z -> pdf(Gumbel(), 40.0) / φ(z)),
        (Normal(0.2, 1.3), StdLogistic(), 40.0, y -> exp(-40.0) / (1 + exp(-40.0))^2 / pdf(Normal(0.2, 1.3), y)),
        (Normal(0.2, 1.3), StdLogistic(), -40.0, y -> exp(-40.0) / (1 + exp(-40.0))^2 / pdf(Normal(0.2, 1.3), y)),
    ]
        f = transport_to(ν, μ)
        @test Zygote.gradient(f, x)[1] ≈ ref(f(x))
    end
end


@testset "Zygote gradients of array products" begin
    # Closed forms: ∂/∂x log N(x | μ, σ) = (μ - x) / σ², ∂/∂μ = (x - μ) / σ²,
    # ∂/∂σ = ((x - μ)² / σ² - 1) / σ.
    μ, σ, x = [0.1, 0.2, 0.3], [1.0, 2.0, 0.5], [0.3, -0.8, 1.1]
    pm = productmeasure(asmeasure.(Normal.(μ, σ)))
    @test Zygote.gradient(x -> logdensityof(pm, x), x)[1] ≈ (μ .- x) ./ σ.^2
    @test Zygote.gradient(μ -> logdensityof(productmeasure(asmeasure.(Normal.(μ, σ))), x), μ)[1] ≈ (x .- μ) ./ σ.^2
    @test Zygote.gradient(σ -> logdensityof(productmeasure(asmeasure.(Normal.(μ, σ))), x), σ)[1] ≈ ((x .- μ).^2 ./ σ.^2 .- 1) ./ σ

    X = [0.3 0.5; -0.8 0.1; 1.1 2.0]
    @test Zygote.gradient(X -> sum(logdensities(pm, X)), X)[1] ≈ (μ .- X) ./ σ.^2
    @test Zygote.gradient(μ -> sum(logdensities(productmeasure(asmeasure.(Normal.(μ, σ))), X)), μ)[1] ≈ vec(sum((X .- μ) ./ σ.^2, dims = 2))

    θ = [0.2, 1.3]
    f_fill = θ -> logdensityof(productmeasure(Fill(asmeasure(Normal(θ[1], θ[2])), 3)), x)
    @test Zygote.gradient(f_fill, θ)[1] ≈ [sum(x .- θ[1]) / θ[2]^2, sum(((x .- θ[1]).^2 ./ θ[2]^2 .- 1) ./ θ[2])]

    # An array of powers:
    Y = [0.3 1.1; -0.8 0.4]
    f_pwr = θ -> logdensityof(productmeasure([asmeasure(Normal(θ[1], 1.0))^2, asmeasure(Normal(θ[2], 2.0))^2]), Y)
    @test Zygote.gradient(f_pwr, [0.1, 0.5])[1] ≈ [sum(Y[:, 1] .- 0.1), sum(Y[:, 2] .- 0.5) / 4]

    # Transports of products built inside the differentiated function,
    # x = μ + σ z:
    z, w = [0.3, -0.8, 1.1], [0.2, -0.5, 0.7]
    f_trafo = (μ, σ, z) -> sum(transport_to(productmeasure(asmeasure.(Normal.(μ, σ))), StdNormal()^3)(z) .* w)
    @test all(Zygote.gradient(f_trafo, μ, σ, z) .≈ (w, z .* w, σ .* w))

    # Likelihoods with parameter-dependent marginals:
    ℒ = likelihoodof(μ -> productmeasure(asmeasure.(Normal.(μ, σ))), x)
    @test Zygote.gradient(μ -> logdensityof(ℒ, μ), μ)[1] ≈ (x .- μ) ./ σ.^2
end
