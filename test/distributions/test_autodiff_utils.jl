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
    # Transports into standard measures currently fail in Zygote: the SArray
    # constructor rrule of StaticArrays returns a thunk for the tuple argument
    # on the standard-variate stream path (src/transport.jl, _transport_via_std).
    for (ν, μ, x, broken) in [
        (StdNormal(), Normal(1, 2), 0.7, true), (Normal(1, 2), StdNormal(), 0.7, false),
        (StdNormal(), Gamma(2.0), 1.3, true), (Gamma(2.0), StdNormal(), 0.4, false),
        (StdLogistic(), TDist(3), 0.8, true), (TDist(3), StdUniform(), 0.3, false),
        (Beta(2, 3), Weibull(0.7, 1.3), 0.9, false), (asmeasure(Logistic()), asmeasure(Normal(0.5)), 0.2, false),
    ]
        f = transport_to(ν, μ)
        ∂f(x) = try Zygote.gradient(f, x)[1] catch; NaN end
        @test ∂f(x) ≈ ForwardDiff.derivative(f, x) broken = broken
    end
end
