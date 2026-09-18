using Test

using MeasureBase.Interface: transport_to, test_transport
using MeasureBase: StdUniform, StdExponential, StdLogistic, StdNormal
using MeasureBase: Dirac, Half, restrict, mbind, productmeasure, pushfwd
using MeasureBase: transport_to_std, transport_from_std, transport_from_std_with_rest
using InverseFunctions: inverse
using MeasureBase: weightedmeasure, mcombine
using StaticArrays: SVector
using Static: static
using LogExpFunctions: logit
using ArraysOfArrays: sliced, flatview, fused
using JLArrays

include("testutils.jl")

@testset "transport_to" begin
    for (f, μ) in [
        (logit, StdUniform())
        (log, StdExponential())
        (exp, StdNormal())
    ]
        test_transport(μ, pushfwd(f, μ))
        test_transport(pushfwd(f, μ), μ)
    end

    for μ0 in [StdUniform(), StdExponential(), StdLogistic(), StdNormal()],
        ν0 in [StdUniform(), StdExponential(), StdLogistic(), StdNormal()]

        @testset "transport_to (variations of) $(nameof(typeof(μ0))) to $(nameof(typeof(ν0)))" begin
            test_transport(ν0, μ0)
            test_transport(2.2 * ν0, 2.2 * μ0)
            test_transport(ν0, μ0^1)
            test_transport(ν0^1, μ0)
            test_transport(ν0^3, μ0^3)
            test_transport(ν0^(2, 3, 2), μ0^(3, 4))
            test_transport(2.2 * ν0^(2, 3, 2), 2.2 * μ0^(3, 4))
            @test_throws ArgumentError transport_to(ν0, μ0)(rand(μ0^12))
            @test_throws ArgumentError transport_to(ν0^3, μ0^3)(rand(μ0^(3, 4)))
        end
    end

    @testset "transfrom from/to Dirac" begin
        μ = Dirac(4.2)
        test_transport(StdExponential()^0, μ)
        test_transport(StdExponential()^(0, 0, 0), μ)
        test_transport(μ, StdExponential()^static(0))
        test_transport(μ, StdExponential()^(static(0), static(0)))
        @test_throws ArgumentError transport_to(StdExponential()^1, μ)
        @test_throws ArgumentError transport_to(μ, StdExponential()^1)
    end

    @testset "transport_to autosel" begin
        @test @inferred(transport_to(StdExponential, StdUniform())) ==
              transport_to(StdExponential(), StdUniform())
        @test @inferred(transport_to(StdExponential, StdUniform()^(2, 3))) ==
              transport_to(StdExponential()^6, StdUniform()^(2, 3))
        @test @inferred(transport_to(StdUniform(), StdExponential)) ==
              transport_to(StdUniform(), StdExponential())
        @test @inferred(transport_to(StdUniform()^(2, 3), StdExponential)) ==
              transport_to(StdUniform()^(2, 3), StdExponential()^6)
    end

    # Tail accuracy of transports between standard measures:

    @testset "transports between standard measures" begin
        stds = (StdUniform(), StdExponential(), StdLogistic(), StdNormal())
        for ν in stds, μ in stds
            f = transport_to(ν, μ)
            for x in [rand(μ) for _ in 1:5]
                @test inverse(f)(f(x)) ≈ x
            end
        end

        # Round trips between the unbounded standard measures keep the tails:
        for z in (-37.0, -20.0, -8.0, -6.0, 6.0, 8.0, 20.0, 37.0)
            for ν in (StdExponential(), StdLogistic())
                y = transport_to(ν, StdNormal())(z)
                @test isfinite(y)
                @test transport_to(StdNormal(), ν)(y) ≈ z rtol = 1e-8
            end
        end
        for l in (-700.0, -40.0, -8.0, 8.0, 40.0, 700.0)
            y = transport_to(StdExponential(), StdLogistic())(l)
            @test isfinite(y) && y >= 0
            @test transport_to(StdLogistic(), StdExponential())(y) ≈ l rtol = 1e-8
        end
        # The lower tail survives a uniform pivot, the upper tail saturates:
        @test transport_to(StdNormal(), StdUniform())(transport_to(StdUniform(), StdNormal())(-37.0)) ≈ -37.0 rtol = 1e-8
    end

    @testset "scalar and static transports" begin
        f = transport_to(StdNormal(), StdUniform())
        @test @inferred(f(0.3)) isa Float64
        @test allocations_of(f, 0.3) == 0
        g = transport_to(StdExponential()^static(3), StdNormal()^static(3))
        xs = SVector(0.1, -0.4, 2.0)
        @test @inferred(g(xs)) isa SVector{3,Float64}
        @test allocations_of(g, xs) == 0
        @test inverse(g)(g(xs)) ≈ xs
        h = transport_to(StdNormal()^3, StdUniform()^3)
        @test h(Float32[0.1, 0.5, 0.9]) isa Vector{Float32}
    end

    @testset "nested powers" begin
        μ = (StdNormal()^2)^3
        x = rand(μ)
        f = transport_to(StdUniform()^6, μ)
        y = f(x)
        @test y isa AbstractVector{<:Real} && length(y) == 6
        x_reco = inverse(f)(y)
        @test all(map(≈, x_reco, x))
        test_transport(StdExponential()^(3, 2), μ)
    end

    @testset "powers of measures without fast DOF" begin
        f_β(a) = StdNormal()^length(a)
        μ = mbind(f_β, StdUniform()^1, vcat)
        P = μ^2
        x = [rand(μ), rand(μ)]
        z = transport_to(StdUniform()^4, P)(x)
        @test z isa AbstractVector{<:Real} && length(z) == 4
        x_reco = transport_to(P, StdUniform()^4)(z)
        @test x_reco isa AbstractVector && all(map(≈, x_reco, x))
    end

    @testset "Half" begin
        μ = Half(StdNormal())
        test_transport(StdUniform(), μ)
        test_transport(StdLogistic(), μ)
        test_transport(μ, StdNormal())
        @test transport_to(StdUniform(), μ)(0.0) ≈ 0
    end

    @testset "array products of mixed standard measures" begin
        src = productmeasure([StdNormal(), StdExponential()])
        trg = productmeasure([StdUniform(), StdLogistic()])
        @test MeasureBase.preferred_stdmeasure(src) === StdNormal
        f = transport_to(trg, src)
        x = [0.5, 1.0]
        y = f(x)
        @test y ≈ [transport_to(StdUniform(), StdNormal())(0.5), transport_to(StdLogistic(), StdExponential())(1.0)]
        @test inverse(f)(y) ≈ x
        @test_throws ArgumentError transport_to_std(MeasureBase.StdMeasure, StdNormal(), 0.5)

        pm = productmeasure(AbstractMeasure[StdNormal(), StdNormal()^2, Dirac(1.0)])
        xm = [0.5, [0.1, 0.2], 1.0]
        z = transport_to(StdUniform()^3, pm)(xm)
        @test z isa AbstractVector{<:Real} && length(z) == 3
        xm_reco = transport_to(pm, StdUniform()^3)(z)
        @test xm_reco[1] ≈ xm[1] && xm_reco[2] ≈ xm[2] && xm_reco[3] == 1.0
    end

    @testset "measures without standard transport" begin
        μ = restrict(x -> x > 0, StdNormal())
        @test_throws ArgumentError transport_to(StdUniform(), μ)(0.5)
        @test_throws ArgumentError transport_to(μ, StdUniform())(0.5)
    end

    @testset "batched transport" begin
        f = transport_to(StdNormal(), StdUniform())
        X = rand(7)
        @test f.(X) ≈ map(f, X)
        @test inverse(f).(f.(X)) ≈ X
        @test eltype(f.(rand(Float32, 5))) == Float32

        g = transport_to(StdExponential()^3, StdNormal()^3)
        Xn = randn(3, 5)
        Yn = g.(sliced(Xn, Val(1)))
        @test Yn isa AbstractVector && length(Yn) == 5
        @test flatview(Yn) ≈ stack(map(g, eachcol(Xn)))
        @test flatview(g.(Xn)) ≈ flatview(Yn)
        @test flatview(inverse(g).(Yn)) ≈ Xn
        Xv = [randn(3) for _ in 1:4]
        @test g.(Xv) == map(g, Xv)

        h = transport_to(StdUniform()^(2, 3), (StdNormal()^2)^3)
        Xh = randn(2, 3, 4)
        @test flatview(h.(Xh)) ≈ stack([h(Xh[:, :, i]) for i in 1:4])
        @test flatview(fused(inverse(h).(h.(Xh)))) ≈ Xh
        Yh = h.(Xh)
        Xh_reco = inverse(h).(Yh)
        @test Xh_reco[2] == inverse(h)(Yh[2])
        @test Xh_reco[2] isa AbstractVector && length(Xh_reco[2]) == 3 && Xh_reco[2][1] isa AbstractVector

        X3 = randn(3, 4, 5)
        Y3 = g.(X3)
        @test size(Y3) == (4, 5) && size(flatview(Y3)) == (3, 4, 5)
        @test Y3[2, 3] ≈ g(X3[:, 2, 3])
        @test f.(SVector(0.3, 0.6, 0.9)) isa SVector{3,Float64}
        @test g.(Xn .+ 0.0) == g.(Xn)

        P = MeasureBase.ProductMeasure([weightedmeasure(log(i), StdNormal()) for i in 1:3])
        p = transport_to(StdUniform()^3, P)
        Xp = randn(3, 6)
        Yp = p.(sliced(Xp, Val(1)))
        @test flatview(Yp) ≈ stack(map(p, eachcol(Xp)))
        @test flatview(inverse(p).(Yp)) ≈ Xp

        mc = mcombine(vcat, StdNormal()^2, StdUniform()^3)
        c = transport_to(StdExponential()^5, mc)
        Xc = vcat(randn(2, 4), rand(3, 4))
        Yc = c.(sliced(Xc, Val(1)))
        @test flatview(Yc) ≈ stack(map(c, eachcol(Xc)))
        @test flatview(inverse(c).(Yc)) ≈ Xc
        cd = transport_to(mcombine(vcat, Dirac(0.5), StdUniform()^2), StdNormal()^2)
        @test flatview(cd.(randn(2, 3)))[1, :] == fill(0.5, 3)

        pf = transport_to(StdUniform(), pushfwd(exp, StdNormal()))
        Xe = exp.(randn(8))
        @test pf.(Xe) ≈ map(pf, Xe)
        @test inverse(pf).(pf.(Xe)) ≈ Xe

        w = transport_to(StdLogistic()^2, weightedmeasure(0.3, StdNormal()^2))
        Xw = randn(2, 5)
        @test flatview(w.(sliced(Xw, Val(1)))) ≈ stack(map(w, eachcol(Xw)))

        JLArrays.allowscalar(false)
        Xj = JLArray(Xn)
        Yj = g.(sliced(Xj, Val(1)))
        @test flatview(Yj) isa JLArray
        @test Array(flatview(Yj)) ≈ flatview(Yn)
        @test Array(flatview(c.(sliced(JLArray(Xc), Val(1))))) ≈ flatview(Yc)
        Pj = MeasureBase.ProductMeasure(JLArray([weightedmeasure(log(i), StdNormal()) for i in 1:3]))
        pj = transport_to(StdUniform()^3, Pj)
        @test Array(flatview(pj.(sliced(JLArray(Xp), Val(1))))) ≈ flatview(Yp)
        Yej = pf.(JLArray(Xe))
        @test Yej isa JLArray && Array(Yej) ≈ pf.(Xe)
    end

    @testset "transport for products" begin
        test_transport(
            StdUniform()^(2, 2),
            productmeasure((StdExponential(), StdLogistic()^3)),
        )
        test_transport(
            productmeasure((StdExponential(), StdLogistic()^3)),
            StdUniform()^(2, 2),
        )

        test_transport(
            StdUniform()^(2, 2),
            productmeasure((a = StdExponential(), b = StdLogistic()^3)),
        )
        test_transport(
            productmeasure((a = StdExponential(), b = StdLogistic()^3)),
            StdUniform()^(2, 2),
        )
    end
end
