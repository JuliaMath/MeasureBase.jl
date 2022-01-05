export SpikeMixture

# TODO: Add `AbstractSuperposition <: AbstractMeasure`, and make SpikeMixture a
# subtype of this 
struct SpikeMixture{M,W,S} <: AbstractMeasure
    m::M   # parent
    w::W   # parent weight
    s::S   # spike weight
end

SpikeMixture(μ,w) = SpikeMixture(μ, w, static(1.0) - w) 

function Pretty.tile(μ::SpikeMixture)
    Pretty.list_layout(Pretty.tile.([μ.m, μ.w]), prefix="SpikeMixture")
end

# TODO: Should this base measure be local? 
@inline function basemeasure(μ::SpikeMixture)
    # Compare formula (1.4) in Joris Bierkens, Sebastiano Grazzi, Frank van der Meulen, Moritz Schauer:
    # Sticky PDMP samplers for sparse and local inference problems. 2020. [https://arxiv.org/abs/2103.08478].
    SpikeMixture(basemeasure(μ.m), static(1.0), static(1.0))
end


@inline function logdensity_def(μ::SpikeMixture, x)
    if iszero(x)
        return log(μ.s) 
    else
        return log(μ.w) + logdensity_def(μ.m, x)
    end
end

function gentype(μ::SpikeMixture)
    gentype(μ.m)
end

function Base.rand(rng::AbstractRNG, T::Type, μ::SpikeMixture)
    return (rand(rng, T) < μ.w) * rand(rng, T, μ.m)
end

testvalue(μ::SpikeMixture) = testvalue(μ.m)
