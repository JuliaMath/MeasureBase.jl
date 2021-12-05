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

function test_interface(μ)
    ###########################################################################
    # basemeasure_depth

    static_depth = basemeasure_depth(μ) 
    dynamic_depth = dynamic_basemeasure_depth(μ)

    if static_depth > dynamic_depth
        @warn "basemeasure_depth($μ) greater than requirement, could add some overhead"
    end
    @test static_depth ≥ dynamic_depth

    ###########################################################################
    # basemeasure_type

    @test (basemeasure_type ∘ typeof)(μ) == (typeof ∘ basemeasure)(μ)
end