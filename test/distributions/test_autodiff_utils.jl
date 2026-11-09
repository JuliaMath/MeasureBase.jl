# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using LinearAlgebra
using Distributions, ArraysOfArrays
using FillArrays: Fill
import ForwardDiff, Zygote
using SpecialFunctions: digamma
using MeasureBase: StdNormal, StdUniform, StdLogistic, transport_to, mbind
using MeasureBase: logdensities, transports, pullbck


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


@testset "ForwardDiff through the SMFs of generic distributions" begin
    # The cdf and quantile of TDist in Distributions don't take dual numbers:
    D = ForwardDiff.derivative
    t = TDist(3.0)
    @test D(x -> MeasureBase.smf(t, x), -2.0) ≈ pdf(t, -2.0)
    @test D(x -> MeasureBase.logcsmf(t, x), 5.0) ≈ -pdf(t, 5.0) / ccdf(t, 5.0)
    @test D(p -> MeasureBase.invsmf(t, p), 0.2) ≈ 1 / pdf(t, quantile(t, 0.2))
    @test D(q -> MeasureBase.invcsmf(t, q), 0.01) ≈ -1 / pdf(t, cquantile(t, 0.01))
    # Half transports via the SMFs of its parent, its density is 2 pdf(t, x):
    φ(z) = exp(-z^2 / 2) / sqrt(2π)
    h = MeasureBase.Half(asmeasure(t))
    z_of_x = transport_to(StdNormal(), h)
    @test D(z_of_x, 0.7) ≈ 2 * pdf(t, 0.7) / φ(z_of_x(0.7))
    x_of_z = transport_to(h, StdNormal())
    @test D(x_of_z, 0.3) ≈ φ(0.3) / (2 * pdf(t, x_of_z(0.3)))
end


@testset "shape derivatives for large shapes" begin
    # ∂P(α, x)/∂α of the regularized lower incomplete gamma function near the
    # median, references from the power series in 160-bit arithmetic:
    for (α, x, ref) in [
        (1e8, 9.997e7, -4.4285249615234345e-7),
        (1e8, 1.0e8, -3.989422807338846e-5),
        (1e8, 1.0003e8, -4.4351727339494425e-7),
        (1e9, 9.99905131670195e8, -1.4011411477882818e-7),
        (1e9, 1.0e9, -1.2615662611152106e-5),
        (1e9, 1.000094868329805e9, -1.401805925048155e-7),
        (1e10, 9.9997e9, -4.431516027665192e-8),
        (1e10, 1.0e10, -3.989422804047572e-6),
        (1e10, 1.00003e10, -4.4321808049267906e-8),
    ]
        @test MeasureBase._gamma_cdf_dshape(α, x) ≈ ref rtol = 1e-10
    end
    # NaN where the series would need more terms than the budget, and for
    # infinite shapes:
    @test isnan(MeasureBase._gamma_cdf_dshape(1e13, 1e13))
    @test isnan(MeasureBase._gamma_cdf_dshape(Inf, 1.0))
    @test all(isnan, MeasureBase._beta_cdf_dshape(Inf, 2.0, 0.5))
    # ∂I_x(α, α)/∂α = -∂I_x(α, α)/∂β at x = 1/2, references from the
    # continued fraction in 160-bit arithmetic:
    for (α, ref) in [(1e4, -0.002821006688412662), (1e6, -0.00028209485054363567), (1e8, -2.8209479236157564e-5)]
        ∂α, ∂β = MeasureBase._beta_cdf_dshape(α, α, 0.5)
        @test ∂α ≈ ref rtol = 1e-9
        @test ∂β ≈ -ref rtol = 1e-9
    end
end


@testset "AD through batched kernels of array-variate products" begin
    # The batched kernels must not accumulate into result buffers, which
    # Zygote can't differentiate and which can't hold dual numbers of
    # captured parameters. Closed forms:
    # ∇ₓ log N(x | μ, I) = μ - x, ∇_μ = x - μ.
    mvn_prod = θ -> productmeasure([asmeasure(MvNormal(θ[1:2], I)), asmeasure(MvNormal(θ[3:4], I))])
    θ = [1.0, 1.0, 0.0, 0.0]
    P = mvn_prod(θ)
    X = reshape(collect(0.1:0.1:1.2), 2, 2, 3)
    μ = reshape(θ, 2, 2)
    @test ForwardDiff.gradient(X -> sum(logdensities(P, X)), X) ≈ μ .- X
    @test Zygote.gradient(X -> sum(logdensities(P, X)), X)[1] ≈ μ .- X
    ∇θ = vec(sum(X .- μ, dims = 3))
    @test ForwardDiff.gradient(θ -> sum(logdensities(mvn_prod(θ), X)), θ) ≈ ∇θ
    @test Zygote.gradient(θ -> sum(logdensities(mvn_prod(θ), X)), θ)[1] ≈ ∇θ
    @test eltype(logdensities(P, Float32.(X))) == typeof(logdensityof(P, Float32.(X[:, :, 1])))

    # Transports from StdNormal^4 shift by μ, so the pullback of P is StdNormal^4:
    f = transport_to(P, StdNormal()^4)
    Z = reshape(collect(0.1:0.1:1.2), 4, 3)
    @test ForwardDiff.gradient(Z -> sum(sum, transports(f, Z)), Z) ≈ ones(4, 3)
    @test Zygote.gradient(Z -> sum(sum, transports(f, Z)), Z)[1] ≈ ones(4, 3)
    @test Zygote.gradient(z -> sum(sum, f(z)), Z[:, 1])[1] ≈ ones(4)
    @test Zygote.gradient(z -> logdensityof(pullbck(f, P), z), Z[:, 1])[1] ≈ -Z[:, 1]

    # Dirichlet marginals: ∂/∂xᵢ = (αᵢ - 1) / xᵢ, ∂/∂αᵢ = log(xᵢ) + ψ(∑α) - ψ(αᵢ):
    dir_prod = α -> productmeasure([asmeasure(Dirichlet(α)), asmeasure(Dirichlet(2 .* α))])
    α = [2.0, 3.0, 4.0]
    XD = cat([0.2 0.3; 0.3 0.3; 0.5 0.4], [0.1 0.6; 0.6 0.2; 0.3 0.2]; dims = 3)
    A = hcat(α, 2 .* α)
    @test ForwardDiff.gradient(X -> sum(logdensities(dir_prod(α), X)), XD) ≈ (A .- 1) ./ XD
    @test Zygote.gradient(X -> sum(logdensities(dir_prod(α), X)), XD)[1] ≈ (A .- 1) ./ XD
    ∇α = vec(sum(log.(XD[:, 1, :]) .+ digamma(sum(α)) .- digamma.(α), dims = 2)) .+
        2 .* vec(sum(log.(XD[:, 2, :]) .+ digamma(2 * sum(α)) .- digamma.(2 .* α), dims = 2))
    @test ForwardDiff.gradient(α -> sum(logdensities(dir_prod(α), XD)), α) ≈ ∇α
end
