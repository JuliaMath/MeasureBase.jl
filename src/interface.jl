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

using JET

function test_interface(μ::M) where {M}
    @eval begin
        μ = $μ
        @testset "$μ" begin
            μ = $μ
            M = $M
           
            ###########################################################################
            # basemeasure_depth
            JET.@test_call basemeasure_depth(μ) 
            static_depth = @inferred basemeasure_depth(μ) 

            dynamic_depth = dynamic_basemeasure_depth(μ)

            if static_depth > dynamic_depth
                @warn "basemeasure_depth($μ) greater than requirement, could add some overhead"
            end
            @test static_depth ≥ dynamic_depth

            ###########################################################################
            # testvalue, logdensityof

            x = @inferred testvalue(μ)
            JET.@test_opt basemeasure(μ, x)
            β = basemeasure(μ, x)

            JET.@test_opt logdensityof(μ, x)
            ℓμ = logdensityof(μ, x)

            JET.@test_opt logdensityof(β, x)
            ℓβ = logdensityof(β, x)

            @test ℓμ ≈ logdensity_def(μ, x) + ℓβ
        end
    end
end
