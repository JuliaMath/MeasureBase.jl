using Test

using MeasureBase.Interface: vartransform, test_vartransform
using MeasureBase: StdUniform, StdExponential, StdLogistic


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
end
