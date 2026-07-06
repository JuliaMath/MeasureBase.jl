# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using Random
import Mooncake
import ForwardDiff

using MeasureBase
using MeasureBase: transport_to
using MeasureBase: isneginf, isposinf, _adignore_call
using MeasureBase: check_dof, require_insupport, _origin_depth
using MeasureBase: logdensityof_rt

_mooncake_gradient(f, x) = Mooncake.value_and_gradient!!(
    Mooncake.prepare_gradient_cache(f, x), f, x
)[2][2]

@testset "Mooncake AD rules" begin
    @test Base.get_extension(MeasureBase, :MeasureBaseMooncakeExt) isa Module

    @testset "zero-derivative primitives" begin
        rng = Random.Xoshiro(789990641)
        Mooncake.TestUtils.test_rule(rng, isneginf, 0.5; is_primitive = true)
        Mooncake.TestUtils.test_rule(rng, isposinf, 0.5; is_primitive = true)
        Mooncake.TestUtils.test_rule(rng, _adignore_call, () -> 42.0; is_primitive = true)
        Mooncake.TestUtils.test_rule(rng, check_dof, StdNormal(), StdUniform(); is_primitive = true)
        Mooncake.TestUtils.test_rule(rng, require_insupport, StdNormal(), 0.5; is_primitive = true)
        Mooncake.TestUtils.test_rule(rng, _origin_depth, StdNormal(); is_primitive = true)
        Mooncake.TestUtils.test_rule(rng, logdensityof_rt, StdNormal(), 0.5; is_primitive = true)
    end

    @testset "@_adignore is ignored" begin
        f_adignore(x) = (MeasureBase.@_adignore x^3; x^2)
        @test _mooncake_gradient(f_adignore, 3.0) ≈ 6.0
    end

    @testset "logdensityof gradients" begin
        x = [0.1, -0.2, 0.3]
        f_ld = x -> logdensityof(StdNormal()^3, x)
        @test _mooncake_gradient(f_ld, x) ≈ ForwardDiff.gradient(f_ld, x)

        f_ldu = x -> logdensityof(StdExponential()^3, x)
        @test _mooncake_gradient(f_ldu, abs.(x)) ≈ ForwardDiff.gradient(f_ldu, abs.(x))
    end

    @testset "transport gradients" begin
        x = [0.1, -0.2, 0.3]
        f_t = x -> sum(transport_to(StdUniform()^3, StdNormal()^3)(x))
        @test _mooncake_gradient(f_t, x) ≈ ForwardDiff.gradient(f_t, x)

        u = [0.3, 0.5, 0.7]
        f_ti = u -> sum(transport_to(StdNormal()^3, StdUniform()^3)(u))
        @test _mooncake_gradient(f_ti, u) ≈ ForwardDiff.gradient(f_ti, u)
    end
end
