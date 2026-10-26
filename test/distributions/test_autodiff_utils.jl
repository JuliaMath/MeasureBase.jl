# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using LinearAlgebra
using Distributions, ArraysOfArrays
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
