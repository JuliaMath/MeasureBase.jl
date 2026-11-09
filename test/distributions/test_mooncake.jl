# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using Random, LinearAlgebra
using Distributions
import Mooncake
import ForwardDiff
using ChainRulesTestUtils: test_rrule

using MeasureBase
using MeasureBase: transport_to, transport_def, asmeasure
using MeasureBase: StdUniform, StdNormal

_mooncake_gradient(f, x) = Mooncake.value_and_gradient!!(
    Mooncake.prepare_gradient_cache(f, x), f, x
)[2][2]

_test_gradient(f, x::Real) = @test _mooncake_gradient(f, x) ≈ ForwardDiff.derivative(f, x)
_test_gradient(f, x::AbstractVector) = @test _mooncake_gradient(f, x) ≈ ForwardDiff.gradient(f, x)

@testset "Mooncake AD with Distributions" begin
    @test Base.get_extension(MeasureBase, :MeasureBaseDistributionsMooncakeExt) isa Module

    @testset "zero-derivative primitives" begin
        rng = Random.Xoshiro(789990641)
        Mooncake.TestUtils.test_rule(
            rng, MeasureBase._dist_params_numtype, Normal(0.2, 1.3);
            is_primitive = true,
        )
    end

    @testset "univariate transport gradients" begin
        _test_gradient(x -> transport_def(StdUniform(), Normal(1.0, 2.0), x), 0.5)
        _test_gradient(u -> transport_def(Normal(1.0, 2.0), StdUniform(), u), 0.3)
        _test_gradient(u -> transport_def(Beta(2.0, 3.0), StdUniform(), u), 0.3)
        _test_gradient(x -> transport_def(StdUniform(), Gamma(2.0, 1.0), x), 0.7)
        _test_gradient(x -> transport_def(StdUniform(), truncated(Normal(0.3, 1.2), -0.5, 1.5), x), 0.4)
        _test_gradient(x -> transport_def(StdUniform(), 2.0 * Weibull(0.7) + 1.0, x), 3.0)
        _test_gradient(x -> transport_def(StdNormal(), Uniform(), x), 0.4)
    end

    @testset "multivariate transport gradients" begin
        mvn = MvNormal([0.3, -2.9], [1.7 0.5; 0.5 2.3])
        _test_gradient(x -> sum(transport_to(StdNormal()^2, mvn)(x)), [0.1, -2.0])
        _test_gradient(y -> sum(transport_to(mvn, StdNormal()^2)(y)), [0.2, 0.7])

        pd = product_distribution([Weibull(0.7), Exponential(1.3), Normal(0.5, 2.0)])
        _test_gradient(x -> sum(transport_to(StdNormal()^3, asmeasure(pd))(x)), [0.4, 0.8, 1.5])

        dirich = Dirichlet([2.0, 3.0, 4.0])
        _test_gradient(u -> MeasureBase.transport_from_std(StdUniform, dirich, u)[1], [0.3, 0.7])
        _test_gradient(x -> sum(MeasureBase.transport_to_std(StdUniform, dirich, vcat(x, 1 - sum(x)))), [0.28, 0.23])
    end

    @testset "incomplete gamma and beta function rules" begin
        # Against finite differences, with respect to all arguments (the
        # rules are reverse-mode rules imported from ChainRules):
        rng = Random.Xoshiro(468128163)
        for (f, args) in [
            (MeasureBase._gamma_cdf, (2.5, 1.3)),
            (MeasureBase._gamma_ccdf, (2.5, 1.3)),
            (MeasureBase._gamma_quantile, (0.7, 0.3)),
            (MeasureBase._gamma_cquantile, (0.7, 0.3)),
            (MeasureBase._beta_cdf, (2.0, 3.0, 0.3)),
            (MeasureBase._beta_quantile, (0.6, 30.0, 0.8)),
        ]
            Mooncake.TestUtils.test_rule(rng, f, args...; is_primitive = true, mode = Mooncake.ReverseMode)
            test_rrule(f, args...)
        end
    end

    @testset "mixed precision" begin
        # Float32 shape parameters with Float64 variates and vice versa,
        # against closed forms: the Gamma(2) density is x e^(-x), the Beta(2, 3)
        # density 12 x (1 - x)^2, and the cdf of Beta(a, 2) is x^a (a + 1 - a x):
        @test _mooncake_gradient(x -> transport_to(StdUniform(), asmeasure(Gamma(2f0, 1f0)))(x), 0.7) ≈ 0.7 * exp(-0.7)
        q = quantile(Gamma(2.0), 0.3)
        @test _mooncake_gradient(u -> transport_to(asmeasure(Gamma(2f0, 1f0)), StdUniform())(u), 0.3) ≈ inv(q * exp(-q)) rtol = 1e-6
        @test _mooncake_gradient(x -> transport_to(StdUniform(), asmeasure(Beta(2f0, 3f0)))(x), 0.3) ≈ 12 * 0.3 * 0.7^2 rtol = 1e-6
        q = quantile(Beta(2.0, 3.0), 0.3)
        @test _mooncake_gradient(u -> transport_to(asmeasure(Beta(2f0, 3f0)), StdUniform())(u), 0.3) ≈ inv(12 * q * (1 - q)^2) rtol = 1e-6

        q = quantile(Beta(2.0, 2.0), 0.3)
        dF_da = q^2 * log(q) * (3 - 2q) + q^2 * (1 - q)
        dF_dq = 6 * q * (1 - q)
        ∂a = _mooncake_gradient(a -> transport_to(asmeasure(Beta(a, 2f0)), StdUniform())(0.3), 2f0)
        @test ∂a isa Float32
        @test ∂a ≈ -dF_da / dF_dq rtol = 1e-5
        P = a -> cdf(Gamma(a), 0.7)
        @test _mooncake_gradient(a -> transport_to(StdUniform(), asmeasure(Gamma(a, 1.0)))(0.7f0), 2.0) ≈
            (P(2.0 + 1e-5) - P(2.0 - 1e-5)) / 2e-5 rtol = 1e-6
    end

    @testset "logdensityof gradients" begin
        for d in [
            Weibull(0.7, 1.3),
            MixtureModel([Normal(-1.0, 1.0), Normal(2.0, 3.0)], [0.3, 0.7]),
            MixtureModel([Normal(-2.0, 1.0), Normal(0.0, 2.0), Normal(3.0, 1.0)], [0.2, 0.5, 0.3]),
        ]
            m = asmeasure(d)
            _test_gradient(x -> logdensityof(m, x[1]), [0.5])
        end

        mvn = MvNormal([0.3, -2.9], [1.7 0.5; 0.5 2.3])
        _test_gradient(x -> logdensityof(asmeasure(mvn), x), [0.1, -2.0])
    end
end
