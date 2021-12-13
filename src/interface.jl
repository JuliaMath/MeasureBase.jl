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
    @testset "$μ" begin
        ###########################################################################
        # basemeasure_type
        
        @test @inferred tbasemeasure_type(M) == @inferred basemeasure_type(μ)

        @info typejoin(basemeasure_type(μ), (typeof ∘ basemeasure)(μ))

        ###########################################################################
        # basemeasure_depth
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

    end
end
