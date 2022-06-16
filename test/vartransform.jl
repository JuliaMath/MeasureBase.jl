using Test

using MeasureBase: vartransform, NoVarTransform
using DensityInterface: logdensityof
using InverseFunctions: inverse
using ChangesOfVariables: with_logabsdet_jacobian

@testset "vartransform" begin
    function test_transform_and_back(ν, μ)
        @testset "vartransform powers from $(nameof(typeof(μ))) to $(ν)" begin
            x = rand(μ)
            @test !(@inferred(vartransform(ν, μ)(x)) isa NoVarTransform)
            f = vartransform(ν, μ)
            y = f(x)
            @test @inferred(inverse(f)(y)) ≈ x
            @test @inferred(with_logabsdet_jacobian(f, x)) isa Tuple{typeof(y),Real}
            @test @inferred(with_logabsdet_jacobian(inverse(f), y)) isa Tuple{typeof(x),Real}
            y2, ladj_fwd = with_logabsdet_jacobian(f, x)
            x2, ladj_inv = with_logabsdet_jacobian(inverse(f), y)
            @test x ≈ x2
            @test y ≈ y2
            @test ladj_fwd ≈ - ladj_inv
            @test ladj_fwd ≈ logdensityof(μ, x) - logdensityof(ν, y)
        end
    end

    for μ in [StdUniform(), StdExponential(), StdLogistic()], ν in [StdUniform(), StdExponential(), StdLogistic()]
        @testset "vartransform powers of $(nameof(typeof(μ))) to $(ν)" begin
            test_transform_and_back(ν, μ)
            test_transform_and_back(ν, μ^1)
            test_transform_and_back(ν^1, μ)
            test_transform_and_back(ν^3, μ^3)
            test_transform_and_back(ν^(2,3,2), μ^(3,4))
            @test_throws ArgumentError vartransform(ν, μ)(rand(μ^12))
            @test_throws ArgumentError vartransform(ν^3, μ^3)(rand(μ^(3,4)))
        end
    end
end
