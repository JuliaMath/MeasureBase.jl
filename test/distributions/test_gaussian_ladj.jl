# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test,
    LinearAlgebra, Distributions, ChangesOfVariables, InverseFunctions, ArraysOfArrays
using MeasureBase
import ForwardDiff
using Adapt: adapt
using JLArrays: JLArray

@testset "Affine Gaussian Jacobians" begin
    for T in (Float32, Float64)
        L = T[2 0; 1 3]
        f = transport_to(asmeasure(MvNormal(T[0.25, -0.5], L * L')), StdNormal()^2)
        X = T[-1 0 2; 2 -0.5 0.3]
        Y, ladjs = transports_with_ladj(f, X)
        @test flatview(Y) == flatview(transports(f, X))
        @test ladjs ≈ fill(T(log(6)), 3) rtol=8eps(T)
        @test last(with_logabsdet_jacobian(inverse(f), first(Y))) ≈ -T(log(6)) rtol=8eps(T)
        _, huge_ladj = with_logabsdet_jacobian(f, fill(sqrt(floatmax(T)), 2))
        @test huge_ladj ≈ T(log(6)) rtol=8eps(T)
    end

    source = asmeasure(MvNormal([0.2, -0.3], [4.0 2; 2 2]))
    target = asmeasure(MvNormal([-0.1, 0.5], [1.0 -1; -1 10]))
    f = transport_to(target, source)
    @test last(with_logabsdet_jacobian(f, [0.3, 0.4])) ≈ log(3/2)
    _, mutable_ladjs = transports_with_ladj(f, zeros(2, 2))
    mutable_ladjs .+= 1
    @test mutable_ladjs ≈ fill(log(3/2) + 1, 2)
    Y, ladjs = transports_with_ladj(f, zeros(2, 0))
    @test isempty(Y) && isempty(ladjs)

    mixed = transport_to(
        asmeasure(MvNormal(zeros(Float32, 2), Float32[4 2; 2 10])),
        StdNormal()^2,
    )
    @test last(with_logabsdet_jacobian(mixed, [0.2, -0.3])) ≈ log(6) rtol=8eps(Float64)

    diagonal = MvNormal(zeros(2), Diagonal([1.0, 4.0]))
    changing = transport_to(asmeasure(diagonal), StdNormal()^2)
    before = last(with_logabsdet_jacobian(changing, zeros(2)))
    diagonal.Σ.diag[2] = 9.0
    @test [before, last(with_logabsdet_jacobian(changing, zeros(2)))] ≈ [log(2), log(3)]
    scalar = transport_to(asmeasure(MvNormal(zeros(2), 2.0I)), StdNormal()^2)
    @test last(with_logabsdet_jacobian(scalar, zeros(2))) ≈ log(2)

    for x in ([Inf, 0.0], [NaN, 0.0])
        y, ladj = with_logabsdet_jacobian(f, x)
        @test ladj == 0
    end
    @test ForwardDiff.derivative(t -> last(with_logabsdet_jacobian(f, [t, 0.2])), 0.4) ≈ 0 atol=64eps(
        Float64,
    )

    X = JLArray([0.3 0.5; 0.4 -0.1])
    _, device_ladjs = transports_with_ladj(adapt(JLArray, f), X)
    @test device_ladjs isa JLArray
    @test Array(device_ladjs) ≈ fill(log(3/2), 2)
end
