module Interface

using Reexport

@reexport using MeasureBase

using MeasureBase: basemeasure_depth, proxy, istrue
using MeasureBase: insupport, basemeasure_sequence, commonbase
using MeasureBase: transport_to, NoTransport

using DensityInterface: logdensityof
using InverseFunctions: inverse
using ChangesOfVariables: with_logabsdet_jacobian
using Tricks: static_hasmethod
using IntervalSets: Interval

export test_interface
export test_transport
export test_smf
export basemeasure_depth
export proxy
export insupport
export basemeasure_sequence
export commonbase

using Test

function dynamic_basemeasure_depth(μ::M) where {M}
    if static_hasmethod(proxy, Tuple{M})
        π = proxy(μ)
        if static_hasmethod(basemeasure, Tuple{typeof(π)})
            basemeasure(π) == basemeasure(μ) && return dynamic_basemeasure_depth(π)
        end
    end
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

            @test AbstractMeasure(μ) isa AbstractMeasure

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

            x = @inferred testvalue(Float64, μ)
            β = @inferred basemeasure(μ, x)

            ℓμ = @inferred logdensityof(μ, x)
            ℓβ = @inferred logdensityof(β, x)

            @test ℓμ ≈ logdensity_def(μ, x) + ℓβ

            @test logdensity_def(μ, testvalue(Float64, μ)) isa Real
        end
    end
end

function test_transport(ν, μ)
    supertype(x) = Any
    supertype(x::Real) = Real
    supertype(x::AbstractArray{<:Real,N}) where {N} = AbstractArray{<:Real,N}

    structisapprox(a, b) = isapprox(a, b)
    function structisapprox(a::NTuple{N,Any}, b::NTuple{N,Any}) where {N}
        all(map(structisapprox, a, b))
    end
    function structisapprox(a::NamedTuple{names}, b::NamedTuple{names}) where {names}
        all(map(structisapprox, values(a), values(b)))
    end

    @testset "transport_to $μ to $ν" begin
        x = rand(μ)
        @test !(@inferred(transport_to(ν, μ)(x)) isa NoTransport)
        f = transport_to(ν, μ)
        y = f(x)
        @test structisapprox(@inferred(inverse(f)(y)), x)
        @test @inferred(with_logabsdet_jacobian(f, x)) isa Tuple{supertype(y),Real}
        @test @inferred(with_logabsdet_jacobian(inverse(f), y)) isa Tuple{supertype(x),Real}
        y2, ladj_fwd = with_logabsdet_jacobian(f, x)
        x2, ladj_inv = with_logabsdet_jacobian(inverse(f), y)
        @test structisapprox(x, x2)
        @test structisapprox(y, y2)
        @test isapprox(ladj_fwd, -ladj_inv, atol = 1e-10)
        @test ladj_fwd ≈ logdensityof(μ, x) - logdensityof(ν, y)
    end
end

function test_smf(μ, n = 100)
    @testset "smf($μ)" begin
        # Get `n` sorted uniforms in O(n) time
        p = rand(n)
        p .+= 0:n-1
        p .*= inv(n)

        F(x) = smf(μ, x)
        Finv(p) = invsmf(μ, p)

        @assert issorted(p)
        x = invsmf.(μ, p)
        @test issorted(x)
        @test all(istrue ∘ insupport(μ), x)

        @test all((Finv ∘ F).(x) .≈ x)

        for j in 1:n
            a = rand()
            b = rand()
            a, b = minmax(a, b)
            x = Finv(a)
            y = Finv(b)
            @test μ(Interval{:open,:closed}(x, y)) ≈ (F(y) - F(x))
        end
    end
end

end # module Interface
