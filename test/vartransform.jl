using Test

using MeasureBase.Interface: vartransform, test_vartransform
using MeasureBase: StdUniform, StdExponential, StdLogistic
using MeasureBase: Dirac


@testset "vartransform" begin
    for μ0 in [StdUniform(), StdExponential(), StdLogistic()], ν0 in [StdUniform(), StdExponential(), StdLogistic()]
        @testset "vartransform (variations of) $(nameof(typeof(μ0))) to $(nameof(typeof(ν0)))" begin
            test_vartransform(ν0, μ0)
            test_vartransform(2.2 * ν0, 3 * μ0)
            test_vartransform(ν0, μ0^1)
            test_vartransform(ν0^1, μ0)
            test_vartransform(ν0^3, μ0^3)
            test_vartransform(ν0^(2,3,2), μ0^(3,4))
            test_vartransform(2.2 * ν0^(2,3,2), 3 * μ0^(3,4))
            @test_throws ArgumentError vartransform(ν0, μ0)(rand(μ0^12))
            @test_throws ArgumentError vartransform(ν0^3, μ0^3)(rand(μ0^(3,4)))
        end
    end

    @testset "transfrom from/to Dirac" begin
        μ = Dirac(4.2)
        test_vartransform(StdExponential()^0, μ)
        test_vartransform(StdExponential()^(0,0,0), μ)
        test_vartransform(μ, StdExponential()^static(0))
        test_vartransform(μ, StdExponential()^(static(0),static(0)))
        @test_throws ArgumentError vartransform(StdExponential()^1, μ)
        @test_throws ArgumentError vartransform(μ, StdExponential()^1)
    end

    @testset "vartransform autosel" begin
        @test @inferred(vartransform(StdExponential, StdUniform())) == vartransform(StdExponential(), StdUniform())
        @test @inferred(vartransform(StdExponential, StdUniform()^(2,3))) == vartransform(StdExponential()^6, StdUniform()^(2,3))
        @test @inferred(vartransform(StdUniform(), StdExponential)) == vartransform(StdUniform(), StdExponential())
        @test @inferred(vartransform(StdUniform()^(2,3), StdExponential)) == vartransform(StdUniform()^(2,3), StdExponential()^6)
    end
end
