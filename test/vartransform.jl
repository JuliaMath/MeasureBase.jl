using Test

using MeasureBase: vartransform, NoVarTransform
using MeasureBase: StdUniform, StdExponential, StdLogistic
using DensityInterface: logdensityof
using InverseFunctions: inverse
using ChangesOfVariables: with_logabsdet_jacobian

@testset "vartransform" begin
    supertype(x::Real) = Real
    supertype(x::AbstractArray{T,N}) where {T,N} = AbstractArray{T,N}

    function test_transform_and_back(ν, μ)
        @testset "vartransform $μ to $ν" begin
            x = rand(μ)
            @test !(@inferred(vartransform(ν, μ)(x)) isa NoVarTransform)
            f = vartransform(ν, μ)
            y = f(x)
            @test @inferred(inverse(f)(y)) ≈ x
            @test @inferred(with_logabsdet_jacobian(f, x)) isa Tuple{supertype(y),Real}
            @test @inferred(with_logabsdet_jacobian(inverse(f), y)) isa Tuple{supertype(x),Real}
            y2, ladj_fwd = with_logabsdet_jacobian(f, x)
            x2, ladj_inv = with_logabsdet_jacobian(inverse(f), y)
            @test x ≈ x2
            @test y ≈ y2
            @test ladj_fwd ≈ - ladj_inv
            @test ladj_fwd ≈ logdensityof(μ, x) - logdensityof(ν, y)
        end
    end

    for μ0 in [StdUniform(), StdExponential(), StdLogistic()], ν0 in [StdUniform(), StdExponential(), StdLogistic()]
        @testset "vartransform (powers of) $(nameof(typeof(μ0))) to $(nameof(typeof(ν0)))" begin
            test_transform_and_back(ν0, μ0)
            test_transform_and_back(ν0, μ0^1)
            test_transform_and_back(ν0^1, μ0)
            test_transform_and_back(ν0^3, μ0^3)
            test_transform_and_back(ν0^(2,3,2), μ0^(3,4))
            @test_throws ArgumentError vartransform(ν0, μ0)(rand(μ0^12))
            @test_throws ArgumentError vartransform(ν0^3, μ0^3)(rand(μ0^(3,4)))
        end
    end
end
