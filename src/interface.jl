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
    ###########################################################################
    # basemeasure_depth
    @testset "$μ" begin
        static_depth = @inferred basemeasure_depth(μ) 

        @test static_depth == tbasemeasure_depth(M)
        dynamic_depth = dynamic_basemeasure_depth(μ)

        if static_depth > dynamic_depth
            @warn "basemeasure_depth($μ) greater than requirement, could add some overhead"
        end
        @test static_depth ≥ dynamic_depth

        ###########################################################################
        # testvalue, logdensityof

        x = @inferred testvalue(μ)
        β = @inferred basemeasure(μ)

        @inferred logdensityof(μ, x)
        @inferred logdensityof(β, x)
        ℓμ = logdensityof(μ, x)
        ℓβ = logdensityof(β, x)

        @test ℓμ ≈ logdensity_def(μ, x) + ℓβ

        ###########################################################################
        # basemeasure_type

        @test @inferred tbasemeasure_type(M) == @inferred basemeasure_type(μ)
        @test basemeasure_type(μ) == (typeof ∘ basemeasure)(μ)
    end
end