using Test

using MeasureBase.Interface: transport_to, test_transport
using MeasureBase: StdUniform, StdExponential, StdLogistic, StdNormal
using MeasureBase: Dirac
using LogExpFunctions: logit

using ChainRulesTestUtils

@testset "transport_to" begin
    test_rrule(MeasureBase._origin_depth, pushfwd(exp, StdUniform()))

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

    @testset "transport for products" begin
        test_transport(StdUniform()^(2, 2), productmeasure((StdExponential(), StdLogistic()^3)))
        test_transport(productmeasure((StdExponential(), StdLogistic()^3)), StdUniform()^(2, 2))

        test_transport(StdUniform()^(2, 2), productmeasure((a = StdExponential(), b = StdLogistic()^3)))
        test_transport(productmeasure((a = StdExponential(), b = StdLogistic()^3)), StdUniform()^(2, 2))
    end
end
