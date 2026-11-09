# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential, productmeasure, insupportof, transport_to
using MeasureBase: weightedmeasure, mintegrate_exp, mcombine
using MeasureBase.InverseFunctions: inverse
using StructArrays: StructArray
import StructArrays
using ArraysOfArrays: flatview
using FillArrays: Fill
using Distributions: Normal, cdf
using MeasureBase: Dirac

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

    @testset "tuple batches with arrays of variates and ragged components" begin
        x2 = [0.1, 0.2]
        P1 = productmeasure((StdNormal()^2, StdNormal()))
        ℓ1 = logdensityof(P1, (x2, 0.3))
        @test logdensities(P1, ([x2, x2, x2], [0.3, 0.3, 0.3])) ≈ fill(ℓ1, 3)
        @test logdensities(productmeasure((a = StdNormal()^2, b = StdNormal())), (a = [x2, x2], b = [0.3, 0.3])) ≈ fill(ℓ1, 2)
        Y = rand(P1^3)
        @test logdensities(P1, StructArrays.components(Y)) ≈ logdensityof.(Ref(P1), Y)
        f = transport_to(StdNormal()^3, P1)
        @test flatview(transports(f, StructArrays.components(Y))) ≈ stack(map(f, Y))
        @test_throws ArgumentError logdensities(P1, ([x2, x2, [0.1]], [0.3, 0.3, 0.3]))

        # Value-dependent variate sizes are evaluated variate by variate:
        b = MeasureBase.mbind(a -> StdNormal()^(a > 0 ? 1 : 2), StdNormal(), vcat)
        Pb = productmeasure((a = StdNormal(), b = b))
        xs = [(a = 0.1, b = [0.3, 0.2]), (a = -0.4, b = [-0.3, 0.2, 0.5]), (a = 1.0, b = [0.6, 0.7])]
        ℓs = logdensityof.(Ref(Pb), xs)
        ℓN(x) = -x^2 / 2 - log(2π) / 2
        @test ℓs ≈ [ℓN(x.a) + sum(ℓN, x.b) for x in xs]
        @test logdensities(Pb, xs) ≈ ℓs
        @test logdensities(Pb, StructArray(xs)) ≈ ℓs
        @test logdensities(Pb, (a = map(x -> x.a, xs), b = map(x -> x.b, xs))) ≈ ℓs
        @test logdensityof(Pb^3, xs) ≈ sum(ℓs)
        @test insupportof(Pb^3, xs) isa Union{Bool,MeasureBase.NoFastInsupport}
        g = transport_to(Pb, Pb)
        @test transports(g, xs) == map(g, xs)
        @test_throws ArgumentError logdensities(productmeasure((a = StdNormal(), b = StdNormal()^2)), [(a = 0.1, b = [0.1, 0.2]), (a = 0.1, b = [0.1])])
    end

    @testset "fill columns" begin
        # Constant columns (`rand` gives them for point masses) mix with
        # dense ones:
        Q = productmeasure((a = StdNormal(), b = StdUniform()))
        Xq = StructArray((a = randn(5), b = Fill(0.25, 5)))
        @test logdensities(Q, Xq) ≈ logdensityof.(Ref(Q), Xq)
        fq = transport_to(StdNormal()^2, Q)
        @test flatview(transports(fq, Xq)) ≈ stack(map(fq, Xq))
        Pd = productmeasure((a = StdNormal(), c = Dirac(2.0)))
        Xd = StructArray((a = randn(4), c = Fill(3.0, 4)))
        @test all(==(-Inf), logdensities(Pd, Xd))
        Xd = rand(Pd^4)
        fd = transport_to(StdUniform()^1, Pd)
        U = transports(fd, Xd)
        @test flatview(U) ≈ reshape(cdf.(Normal(), Xd.a), 1, 4)
        Xr = transports(inverse(fd), U)
        @test Xr.c isa Fill && Xr.a ≈ Xd.a
    end

    @test insupportof(StdUniform()^3, [0.1, 0.5, 0.9]) && !insupportof(StdUniform()^3, [0.1, 1.5, 0.9])
    @test insupportof((StdUniform()^2)^3, rand(2, 3)) && !insupportof((StdUniform()^2)^3, fill(2.0, 2, 3))
    @test insupportof(StdUniform()^3, [0.1, 0.5, 0.9]) isa Bool

    @testset "numeric streams of tuple products" begin
        # Numeric vectors are flat streams, for scalar marginals the vector
        # of the marginal variates:
        T = productmeasure((StdNormal(), StdExponential()))
        @test logdensityof(T, [0.3, 0.7]) == logdensityof(T, (0.3, 0.7))
        for P in (Pt, productmeasure((a = StdNormal(), b = StdExponential()^2)))
            x = rand(P)
            s = vcat(values(x)...)
            ℓ_ref = logdensityof(StdNormal(), s[1]) + logdensityof(StdExponential()^2, s[2:3])
            @test logdensityof(P, s) ≈ logdensityof(P, x) ≈ ℓ_ref
            @test MeasureBase.logdensity_def(P, s) ≈ MeasureBase.logdensity_def(P, x)
            f = transport_to(StdUniform()^3, P)
            @test f(s) ≈ f(x) ≈ vcat(cdf(Normal(), s[1]), 1 .- exp.(-s[2:3]))
            X = rand(P^5)
            S = stack(map(x -> vcat(values(x)...), X))
            @test logdensities(P, S) ≈ logdensityof.(Ref(P), X)
            @test flatview(transports(f, S)) ≈ stack(map(f, X))
            for s_bad in (s[1:2], vcat(s, 0.5))
                @test_throws ArgumentError logdensityof(P, s_bad)
                @test_throws ArgumentError f(s_bad)
            end
            @test_throws ArgumentError logdensities(P, S[1:2, :])
        end
    end

    @testset "nested powers of tuple products with array-valued fields" begin
        p = productmeasure((a = StdNormal(), b = StdNormal()^2))
        m = (p^2)^3
        f = transport_to(m, StdNormal()^18)
        Z = reshape(collect(range(-2, 2, length = 72)), 18, 4)
        ℓN(z) = sum(zi -> -zi^2 / 2 - log(2π) / 2, z)
        X = rand(m^4)
        @test size(X) == (4,) && size(X[1]) == (3,) && size(X[1][1]) == (2,)
        @test logdensityof(m^4, X) ≈ sum(x -> logdensityof(m, x), X)
        Y = transports(f, Z)
        @test size(Y) == (4,)
        # The transport between standard normal and normal products only
        # rearranges the variates:
        @test [logdensityof(m, y) for y in Y] ≈ [ℓN(z) for z in eachcol(Z)]
        @test all(i -> logdensityof(m, Y[i]) ≈ logdensityof(m, f(Z[:, i])), 1:4)
        @test inverse(f).(Y) ≈ collect(eachcol(Z))
        Y2, ladjs = transports_with_ladj(f, Z)
        @test isapprox(ladjs, zeros(4), atol = 1e-12)
        q = productmeasure((a = StdNormal(), c = Dirac([1.0, 2.0])))
        Xq = rand(((q^2)^3)^2)
        @test size(Xq) == (2,) && Xq[1][1][1].c == [1.0, 2.0]
    end

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

    @testset "variate and stream rank of powers of tuple products" begin
        T = productmeasure((StdNormal(), StdExponential()))
        flat(x) = reduce(vcat, [vcat(xi...) for xi in x])
        nflat(x) = reduce(vcat, map(flat, x))

        # Variates are struct arrays over the power dimensions, numeric
        # batches are flat streams:
        @test MeasureBase.mspace_ndims(typeof(T^(2, 3))) == 2
        @test logdensities(T, [0.1 0.2; 0.5 0.7]) ≈ [logdensityof(T, (0.1, 0.5)), logdensityof(T, (0.2, 0.7))]

        for P in (T^3, T^(2, 3), weightedmeasure(0.5, T^(2, 3)))
            n = 2 * length(rand(P))
            X = rand(P^4)
            g = transport_to(P, P)
            Y = g.(X)
            @test length(Y) == 4 && all(size(Y[j]) == size(X[j]) && flat(Y[j]) ≈ flat(X[j]) for j in 1:4)
            f = transport_to(StdUniform()^n, P)
            U = rand(StdUniform()^(n, 4))
            Xu = inverse(f).(U)
            @test all(Xu[j] isa StructArray && size(Xu[j]) == size(X[1]) && flat(Xu[j]) ≈ flat(inverse(f)(U[:, j])) for j in 1:4)
            x = rand(P)
            @test f(flat(x)) ≈ f(x)
            @test_throws ArgumentError f(randn(n - 1))
            d = mintegrate_exp(x -> sum(first, x), P)
            @test logdensities(d, X) ≈ [logdensityof(P, X[j]) + sum(first, X[j]) for j in 1:4]
        end

        # Nested powers of tuple products inside streams:
        P2 = (T^3)^2
        x2 = rand(P2)
        m = mcombine(vcat, productmeasure((StdNormal(), P2)), StdNormal())
        s = vcat(0.3, nflat(x2), -0.2)
        @test logdensityof(m, s) ≈ logdensityof(StdNormal(), 0.3) + logdensityof(P2, x2) + logdensityof(StdNormal(), -0.2)
        h = transport_to(StdUniform()^14, m)
        @test inverse(h)(h(s)) ≈ s
    end
end
