using Test

using MeasureBase.Interface: transport_to, test_transport
using MeasureBase: StdUniform, StdExponential, StdLogistic, StdNormal
using MeasureBase: Dirac, Half, restrict, mbind, productmeasure, pushfwd
using MeasureBase: transport_to_std, transport_from_std, transport_from_std_with_rest
using InverseFunctions: inverse
using StaticArrays: SVector
using Static: static
using LogExpFunctions: logit

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
        @test @allocated(f(0.3)) == 0
        g = transport_to(StdExponential()^static(3), StdNormal()^static(3))
        xs = SVector(0.1, -0.4, 2.0)
        @test @inferred(g(xs)) isa SVector{3,Float64}
        @test @allocated(g(xs)) == 0
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

    @testset "measures without standard transport" begin
        μ = restrict(x -> x > 0, StdNormal())
        @test_throws ArgumentError transport_to(StdUniform(), μ)(0.5)
        @test_throws ArgumentError transport_to(μ, StdUniform())(0.5)
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
