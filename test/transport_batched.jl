# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using MeasureBase
using MeasureBase: StdUniform, StdExponential, StdLogistic, StdNormal, Dirac, GenContext
using MeasureBase: productmeasure, pushfwd, mcombine, weightedmeasure, mbind, PushfwdRootMeasure
using MeasureBase: transport_to_std, transport_from_std, transport_to_std_with_rest
using MeasureBase: batched_transport_to_std, batched_transport_from_std
using MeasureBase: batched_transport_to_std_with_rest, batched_transport_from_std_with_rest
using MeasureBase: batched_rand_impl
using InverseFunctions: inverse
using ArraysOfArrays: sliced, flatview
using StaticArrays: SVector
using Distributions: MvNormal, LogNormal, logpdf
using AffineMaps: Mul, MulAdd
using JLArrays

@testset "batched transport" begin
    stdn_to_u = transport_to(StdUniform(), StdNormal())

    @testset "several variates per stream" begin
        X = vcat(randn(2, 6), rand(3, 6))
        Z, R = batched_transport_to_std_with_rest(StdUniform, StdNormal(), X, (2,))
        @test size(Z) == (2, 6) && size(R) == (3, 6)
        @test Z ≈ stdn_to_u.(X[1:2, :])
        Xb, Rb = batched_transport_from_std_with_rest(StdUniform, StdNormal(), Z, (2,))
        @test Xb ≈ X[1:2, :] && size(Rb) == (0, 6)

        # Powers consume their base with their size as multiplicity:
        Zp, Rp = batched_transport_to_std_with_rest(StdUniform, StdNormal()^2, X, ())
        @test Zp ≈ Z && size(Rp) == (3, 6)
        X2 = vcat(X, X)
        Zp2, Rp2 = batched_transport_to_std_with_rest(StdUniform, StdNormal()^2, X2, (2,))
        @test size(Zp2) == (4, 6) && size(Rp2) == (6, 6)
        @test Zp2 ≈ stdn_to_u.(X2[1:4, :])
        Xp2, _ = batched_transport_from_std_with_rest(StdUniform, StdNormal()^2, Zp2, (2,))
        @test Xp2 ≈ reshape(X2[1:4, :], (2, 2, 6))

        # Combined measures split the rows of each variate by component:
        m = mcombine(vcat, StdNormal()^2, StdUniform()^3)
        Z1 = batched_transport_to_std(StdUniform, m, X)
        Zm, Rm = batched_transport_to_std_with_rest(StdUniform, m, X2, (2,))
        @test size(Rm) == (0, 6) && Zm ≈ vcat(Z1, Z1)
        Xm, _ = batched_transport_from_std_with_rest(StdUniform, m, Zm, (2,))
        @test Xm ≈ reshape(X2, (5, 2, 6))

        # Tuple products consume several variates via their degrees of freedom:
        Pt = productmeasure((StdNormal(), StdExponential()^2))
        Zt = rand(6, 4)
        Xt, Rt = batched_transport_from_std_with_rest(StdUniform, Pt, Zt, (2,))
        @test size(Xt[1]) == (2, 4) && size(Xt[2]) == (2, 2, 4) && size(Rt) == (0, 4)
        for j in 1:4, i in 1:2
            a, b = transport_from_std(StdUniform, Pt, Zt[(3i - 2):(3i), j])
            @test_throws ArgumentError transport_from_std(StdUniform, Pt, Zt[:, j])
            @test Xt[1][i, j] ≈ a && Xt[2][:, i, j] ≈ b
        end

        @test_throws ArgumentError batched_transport_to_std_with_rest(StdUniform, StdNormal()^2, X, (3,))
    end

    @testset "tuple and named tuple products" begin
        Pt = productmeasure((StdNormal(), StdExponential()^2))
        Xt = (randn(4), rand(2, 4))
        Zt = batched_transport_to_std(StdUniform, Pt, Xt)
        @test size(Zt) == (3, 4)
        @test Zt ≈ stack([transport_to_std(StdUniform, Pt, (Xt[1][j], Xt[2][:, j])) for j in 1:4])
        Xr = batched_transport_from_std(StdUniform, Pt, Zt)
        @test Xr[1] ≈ Xt[1] && Xr[2] ≈ Xt[2]
        Pn = productmeasure((a = StdNormal(), b = StdExponential()^2))
        Zn = batched_transport_to_std(StdUniform, Pn, (a = Xt[1], b = Xt[2]))
        @test Zn ≈ Zt
        Xn = batched_transport_from_std(StdUniform, Pn, Zn)
        @test Xn.a ≈ Xt[1] && Xn.b ≈ Xt[2]
        @test batched_transport_to_std(StdUniform, Pt, (Xt[1][1], Xt[2][:, 1])) ≈ Zt[:, 1]
        x1 = batched_transport_from_std(StdUniform, Pt, Zt[:, 1])
        @test x1[1] ≈ Xt[1][1] && x1[2] ≈ Xt[2][:, 1]
        @test_throws ArgumentError batched_transport_from_std(StdUniform, Pt, rand(4, 4))
    end

    @testset "array products of array-variate marginals" begin
        P = productmeasure([weightedmeasure(log(i), StdNormal()^2) for i in 1:3])
        X = randn(2, 3, 5)
        Z = batched_transport_to_std(StdUniform, P, X)
        f = transport_to(StdUniform()^6, P)
        @test size(Z) == (6, 5) && Z ≈ stack([f(X[:, :, j]) for j in 1:5])
        @test batched_transport_from_std(StdUniform, P, Z) ≈ X
        Y = f.(sliced(X, Val(2)))
        @test flatview(Y) ≈ Z
        @test flatview(inverse(f).(Y)) ≈ X
        xs = [randn(2) for _ in 1:3]
        z = transport_to_std(StdUniform, P, xs)
        @test z ≈ f(stack(xs))
        xn = transport_from_std(StdUniform, P, z)
        @test length(xn) == 3 && all(xn[i] ≈ xs[i] for i in 1:3)
        @test_throws ArgumentError batched_transport_to_std(StdUniform, P, randn(2, 2, 5))
        @test_throws ArgumentError batched_transport_from_std(StdUniform, P, rand(5, 5))
    end

    @testset "streams with value-dependent sizes" begin
        f_β(a) = StdNormal()^length(a)
        μb = mbind(f_β, StdUniform()^1, vcat)
        m = mcombine(vcat, μb, StdExponential())
        X = vcat(rand(1, 4), randn(1, 4), rand(1, 4))
        Z = batched_transport_to_std(StdUniform, m, X)
        @test size(Z) == (3, 4)
        @test Z ≈ stack([transport_to_std(StdUniform, m, X[:, j]) for j in 1:4])
        @test batched_transport_from_std(StdUniform, m, Z) ≈ X
        P = μb^2
        x = vcat(rand(1), randn(1), rand(1), randn(1), rand(2))
        z, x_μ, x_rest = transport_to_std_with_rest(StdUniform, P, x)
        @test length(z) == 4 && length(x_μ) == 4 && length(x_rest) == 2
        @test z ≈ transport_to_std(StdUniform, P, [x[1:2], x[3:4]])
    end

    @testset "elementwise pushforwards" begin
        νe = pushfwd(Base.BroadcastFunction(exp), StdNormal()^3)
        @test MeasureBase.mspace_ndims(typeof(νe)) == 1
        Ye = exp.(randn(3, 4))
        @test logdensities(νe, Ye) ≈ [logdensityof(νe, Ye[:, j]) for j in 1:4]
        @test logdensityof(νe, Ye[:, 1]) ≈ sum(logpdf.(LogNormal(), Ye[:, 1]))
        fe = transport_to(StdUniform()^3, νe)
        @test flatview(fe.(sliced(Ye, Val(1)))) ≈ stack(map(fe, eachcol(Ye)))
        @test flatview(inverse(fe).(fe.(sliced(Ye, Val(1))))) ≈ Ye
        @test size(batched_rand_impl(GenContext{Float64}(), νe, (5,))) == (3, 5)
        νr = pushfwd(Base.BroadcastFunction(exp), StdNormal()^3, PushfwdRootMeasure())
        @test logdensities(νr, Ye) ≈ [logdensityof(νr, Ye[:, j]) for j in 1:4]
        JLArrays.allowscalar(false)
        @test Array(logdensities(νe, JLArray(Ye))) ≈ logdensities(νe, Ye)
        @test Array(flatview(fe.(sliced(JLArray(Ye), Val(1))))) ≈ flatview(fe.(sliced(Ye, Val(1))))
    end

    @testset "affine pushforwards" begin
        A = [2.0 0.5; 0.0 1.5]
        b = [1.0, -1.0]
        ν = pushfwd(MulAdd(A, b), StdNormal()^2)
        @test MeasureBase.mspace_ndims(typeof(ν)) == 1
        Y = randn(2, 5)
        @test logdensities(ν, Y) ≈ [logpdf(MvNormal(b, A * A'), Y[:, j]) for j in 1:5]
        @test logdensityof(ν, Y[:, 1]) ≈ logpdf(MvNormal(b, A * A'), Y[:, 1])
        f = transport_to(StdUniform()^2, ν)
        @test flatview(f.(sliced(Y, Val(1)))) ≈ stack(map(f, eachcol(Y)))
        @test flatview(inverse(f).(f.(sliced(Y, Val(1))))) ≈ Y
        @test size(batched_rand_impl(GenContext{Float64}(), ν, (7,))) == (2, 7)
        νs = pushfwd(Mul(2.0), StdNormal())
        @test logdensities(νs, Y[1, :]) ≈ logdensityof.(Ref(νs), Y[1, :])
    end

    @testset "single variates through batched forms" begin
        approx(a::Tuple, b::Tuple) = all(map(approx, a, b))
        approx(a, b) = a ≈ b
        for μ in (
            StdNormal(),
            StdNormal()^3,
            productmeasure((StdNormal(), StdExponential()^2)),
            mcombine(vcat, StdNormal()^2, StdUniform()^3),
            weightedmeasure(0.3, StdNormal()^2),
            pushfwd(Base.BroadcastFunction(exp), StdNormal()^2),
        )
            x = rand(μ)
            z = MeasureBase._as_stdstream(transport_to_std(StdUniform, μ, x))
            zb = batched_transport_to_std(StdUniform, μ, x)
            @test zb isa AbstractVector && zb ≈ z
            @test approx(batched_transport_from_std(StdUniform, μ, z), x)
            @test approx(transport_from_std(StdUniform, μ, MeasureBase._chunk_as_variate(μ, z)), x)
        end
        @test batched_transport_to_std(StdUniform, Dirac(1.0), 1.0) == SVector{0,Bool}()
        @test batched_transport_from_std(StdUniform, Dirac(2.0), SVector{0,Bool}()) == 2.0
        @test_throws ArgumentError transport_to_std(StdUniform, StdNormal()^3, randn(3, 2))
    end
end
