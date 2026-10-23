# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential, productmeasure, insupportof, transport_to
using MeasureBase.InverseFunctions: inverse
using StructArrays: StructArray
using ArraysOfArrays: flatview

@testset "structured batches" begin
    Pt = productmeasure((StdNormal(), StdExponential()^2))
    xs = [rand(Pt) for _ in 1:6]
    ℓ_ref = logdensityof.(Ref(Pt), xs)
    @test logdensities(Pt, xs) ≈ ℓ_ref
    @test logdensities(Pt, StructArray(xs)) ≈ ℓ_ref
    X = rand(Pt^6)
    @test X isa StructArray
    @test logdensities(Pt, X) ≈ logdensityof.(Ref(Pt), X)
    @test logdensityof(Pt^6, X) ≈ sum(logdensityof.(Ref(Pt), X))
    @test logdensityof(Pt^6, collect(X)) ≈ logdensityof(Pt^6, X)
    Xm = rand(Pt^(2, 3))
    @test size(Xm) == (2, 3) && logdensityof(Pt^(2, 3), Xm) ≈ sum(logdensityof.(Ref(Pt), Xm))
    @test_throws ArgumentError logdensityof(Pt^5, X)

    Pn = productmeasure((a = StdNormal(), b = StdExponential()^2))
    Xn = rand(Pn^5)
    @test logdensities(Pn, Xn) ≈ logdensityof.(Ref(Pn), Xn)
    @test logdensities(Pn, collect(Xn)) ≈ logdensityof.(Ref(Pn), Xn)

    f = transport_to(StdUniform()^3, Pt)
    Y = f.(X)
    @test Y isa AbstractVector && length(Y) == 6 && flatview(Y) ≈ stack(map(f, X))
    Xr = inverse(f).(Y)
    @test Xr isa StructArray && all(map((a, b) -> all(map(≈, a, b)), Xr, X))
    h = transport_to(Pn, Pt)
    Yn = h.(X)
    @test Yn isa StructArray && Yn[1] isa NamedTuple{(:a, :b)}
    @test all(Yn[i].a ≈ X[i][1] && Yn[i].b ≈ X[i][2] for i in 1:6)
    @test inverse(h).(Yn) isa StructArray

    @test insupportof(StdUniform()^3, [0.1, 0.5, 0.9]) && !insupportof(StdUniform()^3, [0.1, 1.5, 0.9])
    @test insupportof((StdUniform()^2)^3, rand(2, 3)) && !insupportof((StdUniform()^2)^3, fill(2.0, 2, 3))
    @test insupportof(StdUniform()^3, [0.1, 0.5, 0.9]) isa Bool

    @testset "transports of powers of tuple products" begin
        flat(x) = reduce(vcat, [vcat(xi...) for xi in x])
        for P in (productmeasure((StdNormal(), StdExponential()))^3, productmeasure((a = StdNormal(), b = StdExponential()))^3)
            f = transport_to(StdUniform()^6, P)
            finv = inverse(f)
            x = rand(P)
            u = f(x)
            @test length(u) == 6 && all(0 .< u .< 1)
            @test flat(finv(u)) ≈ flat(x)
            U = rand(StdUniform()^(6, 4))
            Y = finv.(U)
            @test length(Y) == 4
            @test all(flat(Y[j]) ≈ flat(finv(U[:, j])) for j in 1:4)
            @test flatview(f.(Y)) ≈ U
        end
    end

    @testset "densities of powers of tuple products" begin
        T = productmeasure((StdNormal(), StdExponential()))
        x = rand(T^3)
        xs = reduce(vcat, [[a, b] for (a, b) in x])
        @test logdensityof(T^3, xs) ≈ logdensityof(T^3, x)
        @test MeasureBase.logdensity_def(T^3, xs) ≈ MeasureBase.logdensity_def(T^3, x)
        @test logdensities(T^3, (randn(3, 2), rand(3, 2))) isa AbstractVector
        @test_throws ArgumentError logdensities(T^3, (randn(5), rand(5)))
        @test_throws ArgumentError logdensities(T^3, (randn(2, 4), rand(2, 4)))
    end
end
