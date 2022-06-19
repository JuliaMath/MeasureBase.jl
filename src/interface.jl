module Interface

using Reexport

@reexport using MeasureBase

using MeasureBase: basemeasure_depth, proxy
using MeasureBase: insupport, basemeasure_sequence, commonbase
using MeasureBase: transport_to, NoVarTransform

using DensityInterface: logdensityof
using InverseFunctions: inverse
using ChangesOfVariables: with_logabsdet_jacobian

export test_interface
export test_vartransform
export basemeasure_depth
export proxy
export insupport
export basemeasure_sequence
export commonbase

using Test

function dynamic_basemeasure_depth(μ)
    β = basemeasure(μ)
    depth = 0
    while μ ≠ β
        (μ, β) = (β, basemeasure(β))
        depth += 1
    end
    return depth
end

function test_interface(μ::M) where {M}
    @eval begin
        μ = $μ
        @testset "$μ" begin
            μ = $μ

            ###########################################################################
            # basemeasure_depth
            static_depth = @inferred basemeasure_depth(μ)

            dynamic_depth = dynamic_basemeasure_depth(μ)

            if static_depth > dynamic_depth
                @warn "basemeasure_depth($μ) greater than requirement, could add some overhead"
            end
            @test static_depth ≥ dynamic_depth

            ###########################################################################
            # testvalue, logdensityof

            x = @inferred testvalue(μ)
            β = @inferred basemeasure(μ, x)

            ℓμ = @inferred logdensityof(μ, x)
            ℓβ = @inferred logdensityof(β, x)

            @test ℓμ ≈ logdensity_def(μ, x) + ℓβ

            @test logdensity_def(μ, testvalue(μ)) isa Real
        end
    end
end


function test_vartransform(ν, μ)
    supertype(x::Real) = Real
    supertype(x::AbstractArray{<:Real,N}) where N = AbstractArray{<:Real,N}

    @testset "transport_to $μ to $ν" begin
        x = rand(μ)
        @test !(@inferred(transport_to(ν, μ)(x)) isa NoVarTransform)
        f = transport_to(ν, μ)
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

end # module Interface
