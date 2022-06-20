using Test

using MeasureBase.Interface: transport_to, test_transport
using MeasureBase: StdUniform, StdExponential, StdLogistic
using MeasureBase: Dirac


@testset "transport_to" begin
    for μ0 in [StdUniform(), StdExponential(), StdLogistic()], ν0 in [StdUniform(), StdExponential(), StdLogistic()]
        @testset "transport_to (variations of) $(nameof(typeof(μ0))) to $(nameof(typeof(ν0)))" begin
            test_transport(ν0, μ0)
            test_transport(2.2 * ν0, 3 * μ0)
            test_transport(ν0, μ0^1)
            test_transport(ν0^1, μ0)
            test_transport(ν0^3, μ0^3)
            test_transport(ν0^(2,3,2), μ0^(3,4))
            test_transport(2.2 * ν0^(2,3,2), 3 * μ0^(3,4))
            @test_throws ArgumentError transport_to(ν0, μ0)(rand(μ0^12))
            @test_throws ArgumentError transport_to(ν0^3, μ0^3)(rand(μ0^(3,4)))
        end
    end

    @testset "transfrom from/to Dirac" begin
        μ = Dirac(4.2)
        test_transport(StdExponential()^0, μ)
        test_transport(StdExponential()^(0,0,0), μ)
        test_transport(μ, StdExponential()^static(0))
        test_transport(μ, StdExponential()^(static(0),static(0)))
        @test_throws ArgumentError transport_to(StdExponential()^1, μ)
        @test_throws ArgumentError transport_to(μ, StdExponential()^1)
    end

    @testset "transport_to autosel" begin
        @test @inferred(transport_to(StdExponential, StdUniform())) == transport_to(StdExponential(), StdUniform())
        @test @inferred(transport_to(StdExponential, StdUniform()^(2,3))) == transport_to(StdExponential()^6, StdUniform()^(2,3))
        @test @inferred(transport_to(StdUniform(), StdExponential)) == transport_to(StdUniform(), StdExponential())
        @test @inferred(transport_to(StdUniform()^(2,3), StdExponential)) == transport_to(StdUniform()^(2,3), StdExponential()^6)
    end
end
