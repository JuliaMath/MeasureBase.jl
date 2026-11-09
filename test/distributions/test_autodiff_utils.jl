# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using LinearAlgebra
using Distributions, ArraysOfArrays
using FillArrays: Fill
import ForwardDiff, Zygote
using MeasureBase: StdNormal, StdUniform, StdLogistic, transport_to, mbind


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


@testset "Zygote gradients of distributions built inside" begin
    # With ForwardDiffPullbacks loaded (as in this test set), the families
    # with native rules and the generic cdf/quantile path. Derivatives of
    # z = Φ⁻¹(F(x)): dz/dx = f(x) / φ(z), dz/dθ = ∂F(x)/∂θ / φ(z).
    φ(z) = exp(-z^2 / 2) / sqrt(2π)
    @test Zygote.gradient(x -> transport_to(StdNormal(), asmeasure(Normal(0.2, 1.3)))(x), 0.7)[1] ≈ 1 / 1.3
    @test Zygote.gradient(a -> transport_to(Normal(a, 1.5), StdUniform())(0.3), 0.5)[1] ≈ 1
    @test Zygote.gradient(θ -> transport_to(StdNormal(), asmeasure(Normal(θ[1], θ[2])))(0.7), [0.2, 1.3])[1] ≈ [-1, -0.5 / 1.3] ./ 1.3

    t_inside = x -> transport_to(StdNormal(), asmeasure(TDist(3.0)))(x)
    t_const = transport_to(StdNormal(), TDist(3.0))
    for f in (t_inside, t_const)
        @test Zygote.gradient(f, 0.7)[1] ≈ pdf(TDist(3.0), 0.7) / φ(f(0.7))
    end
    t_inv = z -> transport_to(asmeasure(TDist(3.0)), StdNormal())(z)
    @test Zygote.gradient(t_inv, 0.7)[1] ≈ φ(0.7) / pdf(TDist(3.0), t_inv(0.7))
    @test Zygote.gradient(x -> transport_to(StdNormal(), asmeasure(Chisq(3.0)))(x), 0.7)[1] ≈
        pdf(Chisq(3.0), 0.7) / φ(transport_to(StdNormal(), Chisq(3.0))(0.7))
    g = a -> transport_to(StdNormal(), asmeasure(Gumbel(a, 2.0)))(0.4)
    @test Zygote.gradient(g, 0.3)[1] ≈ -pdf(Gumbel(0.3, 2.0), 0.4) / φ(g(0.3))

    # Binds with univariate kernels, (a, x) = (z₁, F⁻¹(Φ(z₂))):
    w(y) = y[1] + 2 * y[2]
    bn = z -> w(transport_to(mbind(a -> asmeasure(Normal(a, 1.3)), StdNormal(), tuple), StdNormal()^2)(z))
    @test Zygote.gradient(bn, [0.3, -0.4])[1] ≈ [3.0, 2.6]
    bt = z -> w(transport_to(mbind(a -> asmeasure(TDist(3.0)), StdNormal(), tuple), StdNormal()^2)(z))
    x_t = transport_to(TDist(3.0), StdNormal())(-0.4)
    @test Zygote.gradient(bt, [0.3, -0.4])[1] ≈ [1.0, 2 * φ(-0.4) / pdf(TDist(3.0), x_t)]
end


@testset "nested ForwardDiff through generic transports" begin
    # Second derivatives of x(z) = F⁻¹(Φ(z)) and z(x) = Φ⁻¹(F(x)) from the
    # density f and its log-derivative g = f'/f:
    φ(z) = exp(-z^2 / 2) / sqrt(2π)
    d2x_dz2(f, g, z, x) = (r = φ(z) / f(x); r * (-z - r * g(x)))
    d2z_dx2(f, g, x, z) = (r = f(x) / φ(z); r * (g(x) + z * r))
    D = ForwardDiff.derivative

    t = TDist(3.0)
    f_t, g_t = x -> pdf(t, x), x -> -4x / (3 + x^2)
    x_of_z = transport_to(asmeasure(t), StdNormal())
    z_of_x = transport_to(StdNormal(), asmeasure(t))
    @test D(z -> D(x_of_z, z), 0.3) ≈ d2x_dz2(f_t, g_t, 0.3, x_of_z(0.3))
    @test D(x -> D(z_of_x, x), 0.3) ≈ d2z_dx2(f_t, g_t, 0.3, z_of_x(0.3))

    c = Chisq(3.0)
    z_of_x_c = transport_to(StdNormal(), asmeasure(c))
    @test D(x -> D(z_of_x_c, x), 0.7) ≈ d2z_dx2(x -> pdf(c, x), x -> 1 / (2x) - 1 / 2, 0.7, z_of_x_c(0.7))

    # sum(abs2, x(z)) has the diagonal Hessian 2 (x'² + x x''):
    z = [0.3, -0.2]
    x_of_z2 = transport_to(asmeasure(t)^2, StdNormal()^2)
    H = ForwardDiff.hessian(z -> sum(abs2, x_of_z2(z)), z)
    x = x_of_z2(z)
    @test H ≈ Diagonal(2 .* ((φ.(z) ./ f_t.(x)).^2 .+ x .* d2x_dz2.(f_t, g_t, z, x)))
end
