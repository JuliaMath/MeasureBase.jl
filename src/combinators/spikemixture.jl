export SpikeMixture

# TODO: Add `AbstractSuperposition <: AbstractMeasure`, and make SpikeMixture a
# subtype of this 
struct SpikeMixture{T,S} <: AbstractMeasure
    m::T # parent
    w::S # relative weight of parent
    s::S # scale
end
SpikeMixture(m, w) = SpikeMixture(m, w, one(w))

function Base.show(io::IO, μ::SpikeMixture)
    io = IOContext(io, :compact => true)
    print(io, "(", μ.s * μ.w, "*", string(μ.m), " + ", μ.s * (1 - μ.w), "Dirac(0))")
end

@inline function basemeasure(μ::SpikeMixture)
    # Compare formula (1.4) in Joris Bierkens, Sebastiano Grazzi, Frank van der Meulen, Moritz Schauer:
    # Sticky PDMP samplers for sparse and local inference problems. 2020. [https://arxiv.org/abs/2103.08478].
    ki = (1 / μ.w - 1) / density(μ.m, 0)
    SpikeMixture(basemeasure(μ.m), 1 / (1 + ki), μ.s * (1 + ki))
end

basemeasure_depth(::Type{SpikeMixture{T,S}}) where {T,S} = static(1) + basemeasure_depth(T)

@inline function logdensity(μ::SpikeMixture, x)
    return log(μ.w) + logdensity(μ.m, x)
end

function sampletype(μ::SpikeMixture)
    sampletype(μ.m)
end

function Base.rand(rng::AbstractRNG, T::Type, μ::SpikeMixture)
    μ.s != 1 && throw(ArgumentError("Not a probability measure"))
    return (rand(rng, T) < μ.w) * rand(rng, T, μ.m)
end

testvalue(μ::SpikeMixture) = testvalue(μ.m)
