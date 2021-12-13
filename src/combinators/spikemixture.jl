export SpikeMixture

# TODO: Add `AbstractSuperposition <: AbstractMeasure`, and make SpikeMixture a
# subtype of this 
struct SpikeMixture{T,S} <: AbstractMeasure
    m::T # parent
    w::S # relative weight of parent
end

function Pretty.tile(μ::SpikeMixture)
    result = Pretty.tile(μ.w)
    result *= Pretty.tile(μ.m)
    result *= Pretty.literal(" + ")
    result *= Pretty.tile(1 - μ.w)
    result *= Pretty.literal("Dirac(0)")
    return result
end

# TODO: Should this base measure be local? 
@inline function basemeasure(μ::SpikeMixture)
    # Compare formula (1.4) in Joris Bierkens, Sebastiano Grazzi, Frank van der Meulen, Moritz Schauer:
    # Sticky PDMP samplers for sparse and local inference problems. 2020. [https://arxiv.org/abs/2103.08478].
    SpikeMixture(basemeasure(μ.m), μ.w)
end

tbasemeasure_depth(::Type{SpikeMixture{T,S}}) where {T,S} = static(1) + tbasemeasure_depth(T)

# basemeasure_type(::Type{SpikeMixture{T,S}}) where {T,S} = SpikeMixture{}

function tbasemeasure_type(::Type{SpikeMixture{M, T}}) where {M,T}
    B = tbasemeasure_type(M)
    SpikeMixture{B, T}
end

@inline function logdensity_def(μ::SpikeMixture, x)
    if iszero(x)
        return log1p(-μ.w) 
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
