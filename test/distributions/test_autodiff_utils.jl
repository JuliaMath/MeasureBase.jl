# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using LinearAlgebra
using Distributions, ArraysOfArrays
import ForwardDiff, Zygote


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
