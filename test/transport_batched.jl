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
using ArraysOfArrays: sliced, flatview, nestedview, ArrayOfSimilarArrays
using StaticArrays: SVector, SMatrix, @SMatrix
using MeasureBase: logdensities
using Distributions: MvNormal, LogNormal, Normal, logpdf, cdf
using ChangesOfVariables: with_logabsdet_jacobian
using StructArrays: StructArray, StructArrays
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

    @testset "transports" begin
        f = transport_to(StdUniform()^3, StdNormal()^3)
        X = randn(3, 5)
        Y = transports(f, X)
        @test length(Y) == 5 && flatview(Y) ≈ cdf.(Normal(), X)
        @test transports(f, nestedview(X)) == Y && f.(X) == Y
        @test transports(f, X[:, 1]) ≈ f(X[:, 1])
        g = transport_to(StdUniform(), StdNormal())
        @test transports(g, 0.3) == g(0.3) && transports(g, X) ≈ cdf.(Normal(), X)

        # Results are materialized, not views of temporary arrays:
        fw = transport_to(StdLogistic(), weightedmeasure(0.3, StdNormal()))
        x = randn(5)
        @test transports(fw, x) isa Vector{Float64} && transports(fw, x) ≈ map(fw, x)
        ft = transport_to(productmeasure((StdNormal(), StdExponential())), StdUniform()^2)
        U = rand(2, 5)
        Yt = transports(ft, U)
        @test Yt isa StructArray && all(c -> c isa Array, StructArrays.components(Yt))
        @test all(Yt[j][1] ≈ ft(U[:, j])[1] && Yt[j][2] ≈ ft(U[:, j])[2] for j in 1:5)

        # Fused arguments are materialized once, nested ones into flat storage:
        v = [0.5, 1.0, 1.5]
        @test flatview(f.(nestedview(X) .+ Ref(v))) ≈ cdf.(Normal(), X .+ v)
        @test flatview(f.(X .+ 1)) ≈ cdf.(Normal(), X .+ 1)

        # Tuples and named tuples of batches are batches of tuple variates:
        Pt = productmeasure((StdNormal(), StdExponential()^2))
        Xt = (randn(4), rand(2, 4))
        ht = transport_to(StdUniform()^3, Pt)
        @test flatview(transports(ht, Xt)) ≈ stack([ht((Xt[1][j], Xt[2][:, j])) for j in 1:4])
        @test transports(ht, (Xt[1][1], Xt[2][:, 1])) ≈ ht((Xt[1][1], Xt[2][:, 1]))
        hn = transport_to(StdUniform()^3, productmeasure((a = StdNormal(), b = StdExponential()^2)))
        @test transports(hn, (a = Xt[1], b = Xt[2])) ≈ transports(ht, Xt)

        # Batches that don't fit the variate rank or size throw:
        @test_throws ArgumentError transports(f, randn(2, 5))
        @test_throws ArgumentError transports(f, nestedview(randn(3, 2, 4), 2))
    end

    @testset "results share no memory with the input" begin
        storage(A::AbstractArray{<:Number}) = (A,)
        storage(A::ArrayOfSimilarArrays) = (flatview(A),)
        storage(A::StructArray) = Tuple(Iterators.flatten(map(storage, values(StructArrays.components(A)))))
        storage(X::NamedTuple) = Tuple(Iterators.flatten(map(storage, values(X))))
        P = productmeasure((a = StdNormal(), b = StdLogistic()^2))
        mc = mcombine(vcat, StdNormal()^2, StdNormal()^1)
        for (f, X) in (
            (transport_to(StdNormal(), StdNormal()), randn(5)),
            (transport_to(StdNormal()^3, StdNormal()^3), randn(3, 5)),
            (transport_to(StdNormal()^3, StdNormal()^3), nestedview(randn(3, 5))),
            (transport_to(StdNormal()^3, StdNormal()^3), randn(3)),
            (transport_to(StdLogistic()^(2, 3), StdLogistic()^(2, 3)), randn(2, 3, 4)),
            (transport_to(StdNormal()^2, productmeasure([StdNormal(), StdNormal()])), randn(2, 4)),
            (transport_to(StdNormal()^3, mc), randn(3, 4)),
            (transport_to(mc, mc), randn(3, 4)),
            (transport_to(P, P), StructArray((a = randn(4), b = nestedview(randn(2, 4))))),
            (transport_to(P, P), (a = randn(4), b = randn(2, 4))),
        )
            Y = transports(f, X)
            @test Y !== X
            Y_ref = deepcopy(Y)
            foreach(A -> fill!(A, NaN), storage(X))
            @test isequal(Y, Y_ref)
        end
        f = transport_to(StdNormal()^3, StdNormal()^3)
        JLArrays.allowscalar(false)
        Xj = JLArray(randn(3, 5))
        Yj = transports(f, Xj)
        Yj_ref = Array(flatview(Yj))
        fill!(Xj, NaN)
        @test Array(flatview(Yj)) == Yj_ref
        # Single array variates as well:
        x = randn(3)
        y = f(x)
        @test y !== x && f.(x) !== x
        y_ref = copy(y)
        x .= NaN
        @test y == y_ref
    end

    @testset "batched log-abs-det Jacobians" begin
        # Against the point transports and against the density ratio (the
        # densities of the standard uniform targets are one):
        function test_ladjs(f, X, xs, ladjs_ref)
            Y, ladjs = transports_with_ladj(f, X)
            @test isequal(Y, transports(f, X))
            @test size(ladjs) == size(ladjs_ref) && ladjs ≈ ladjs_ref
            @test isapprox(vec(collect(ladjs)), [last(with_logabsdet_jacobian(f, x)) for x in xs], nans = true)
        end
        normal_logd(x) = sum(logpdf.(Normal(), x))
        X = randn(3, 5)
        test_ladjs(transport_to(StdUniform(), StdNormal()), X[1, :], X[1, :], logpdf.(Normal(), X[1, :]))
        f = transport_to(StdUniform()^3, StdNormal()^3)
        ref = map(normal_logd, eachcol(X))
        test_ladjs(f, X, eachcol(X), ref)
        test_ladjs(f, nestedview(X), eachcol(X), ref)
        @test last(transports_with_ladj(f, X[:, 1])) ≈ normal_logd(X[:, 1])
        @test last(transports_with_ladj(transport_to(StdUniform(), StdNormal()), 0.3)) ≈ logpdf(Normal(), 0.3)
        test_ladjs(transport_to(StdUniform()^(3, 5), StdNormal()^(3, 5)), X, (X,), normal_logd(X))
        # Weights don't change the Jacobians:
        test_ladjs(transport_to(StdUniform(), weightedmeasure(0.7, StdNormal())), X[1, :], X[1, :], logpdf.(Normal(), X[1, :]))
        test_ladjs(transport_to(weightedmeasure(0.7, StdUniform()^3), 2.5 * StdNormal()^3), X, eachcol(X), ref)
        # Tuple products, as tuples of batches and as struct arrays:
        P = productmeasure((a = StdNormal(), b = StdExponential()^2))
        Xt = (a = X[1, :], b = abs.(X[2:3, :]))
        xs = [(a = Xt.a[j], b = Xt.b[:, j]) for j in 1:5]
        ref_t = [normal_logd(x.a) - sum(x.b) for x in xs]
        test_ladjs(transport_to(StdUniform()^3, P), Xt, xs, ref_t)
        test_ladjs(transport_to(StdUniform()^3, P), StructArray(xs), xs, ref_t)
        # Measures containing binds:
        b = mbind(a -> StdNormal()^2, StdNormal(), vcat)
        test_ladjs(transport_to(StdUniform()^3, b), X, eachcol(X), ref)
        test_ladjs(transport_to(StdUniform()^3, b), collect(eachcol(X)), eachcol(X), ref)
        U = rand(4, 5)
        hc = transport_to(mcombine(vcat, b, StdNormal()^1), StdUniform()^4)
        test_ladjs(hc, U, eachcol(U), [-normal_logd(hc(u)) for u in eachcol(U)])
        ht = transport_to(productmeasure((b, StdNormal())), StdUniform()^4)
        test_ladjs(ht, U, eachcol(U), [-normal_logd(y[1]) - normal_logd(y[2]) for y in map(ht, eachcol(U))])
        # Variates outside the support of the source:
        g = transport_to(StdUniform(), StdExponential())
        xe = [-1.0, 0.5, 2.0, NaN]
        @test isequal(last(transports_with_ladj(g, xe)), [last(with_logabsdet_jacobian(g, x)) for x in xe])
        # Static batches stay static, Float32 stays Float32:
        Xs = SMatrix{3,5}(X)
        Ys, ladjs_s = @inferred transports_with_ladj(f, Xs)
        @test flatview(Ys) isa SMatrix{3,5} && ladjs_s isa SVector{5} && ladjs_s ≈ ref
        _, ladjs_32 = transports_with_ladj(f, Float32.(X))
        @test eltype(ladjs_32) == Float32 && isapprox(ladjs_32, ref, rtol = 1e-5)
        # Device arrays:
        JLArrays.allowscalar(false)
        Yj, ladjs_j = transports_with_ladj(f, JLArray(X))
        @test ladjs_j isa JLArray && Array(ladjs_j) ≈ ref && Array(flatview(Yj)) ≈ flatview(transports(f, X))
    end

    @testset "transported batches have the layout of rand" begin
        for ν in (
            productmeasure((StdNormal(), StdExponential())),
            productmeasure((a = StdNormal(), b = StdExponential()^2, c = Dirac(2.0))),
            productmeasure((a = weightedmeasure(0.3, StdNormal()), b = productmeasure((x = StdUniform(),)))),
            productmeasure((StdNormal(), StdExponential()))^3,
            weightedmeasure(0.3, StdNormal()^2),
        )
            n = Int(MeasureBase.getdof(ν))
            Y = transports(transport_to(ν, StdUniform()^n), rand(n, 5))
            @test typeof(Y) == typeof(rand(ν^5))
        end
    end

    @testset "batch layouts" begin
        f = transport_to(StdNormal()^3, StdUniform()^3)
        @test_throws ArgumentError f.(rand(6, 2))
        g = transport_to(StdNormal()^(2, 3), StdUniform()^(2, 3))
        @test_throws ArgumentError g.(rand(3, 2, 4))
        Y = g.(fill(0.5, 2, 3, 4))
        @test size(Y) == (4,) && all(y -> size(y) == (2, 3) && all(iszero, y), Y)

        # Static batches stay static, also nested ones:
        Ys = f.(@SMatrix fill(0.5, 3, 2))
        @test flatview(Ys) isa SMatrix{3,2} && all(iszero, flatview(Ys))
        Xn = SVector(SVector(0.1, 0.2, 0.3), SVector(0.4, 0.5, 0.6))
        @test f.(Xn) ≈ map(f, Xn)
        @test logdensities(StdNormal()^3, Xn) ≈ logdensityof.(Ref(StdNormal()^3), Xn)

        # Flat batches for targets without a variate rank:
        h = transport_to(mbind(x -> StdNormal()^2, StdNormal()), StdNormal()^3)
        X = randn(3, 4)
        @test h.(X) == map(h, eachcol(X))
        @test h.(X[:, 1]) == h(X[:, 1])
    end
end
