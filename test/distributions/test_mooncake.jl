# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using Random, LinearAlgebra
using Distributions
import Mooncake
import ForwardDiff

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
        _test_gradient(x -> transport_def(StdNormal(), StandardDist{Uniform}(), x), 0.4)
    end

    @testset "multivariate transport gradients" begin
        mvn = MvNormal([0.3, -2.9], [1.7 0.5; 0.5 2.3])
        _test_gradient(x -> sum(transport_to(StdNormal()^2, mvn)(x)), [0.1, -2.0])
        _test_gradient(y -> sum(transport_to(mvn, StdNormal()^2)(y)), [0.2, 0.7])

        pd = product_distribution([Weibull(0.7), Exponential(1.3), Normal(0.5, 2.0)])
        _test_gradient(x -> sum(transport_to(StdNormal()^3, asmeasure(pd))(x)), [0.4, 0.8, 1.5])

        dirich = Dirichlet([2.0, 3.0, 4.0])
        _test_gradient(u -> MeasureBase.from_origin(dirich, u)[1], [0.3, 0.7])
        _test_gradient(x -> sum(MeasureBase.to_origin(dirich, vcat(x, 1 - sum(x)))), [0.28, 0.23])
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
