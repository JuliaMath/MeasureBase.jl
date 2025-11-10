export SpikeMixture

# TODO: Add `AbstractSuperposition <: AbstractMeasure`, and make SpikeMixture a
# subtype of this 
struct SpikeMixture{M,W,S} <: AbstractMeasure
    m::M   # parent
    w::W   # parent weight
    s::S   # spike weight
end

SpikeMixture(μ, w) = SpikeMixture(μ, w, static(1.0) - w)

function Pretty.tile(μ::SpikeMixture)
    Pretty.list_layout(Pretty.tile.([μ.m, μ.w]), prefix = "SpikeMixture")
end

# TODO: Should this base measure be local? 
@inline function basemeasure(μ::SpikeMixture)
    # Compare formula (1.4) in Joris Bierkens, Sebastiano Grazzi, Frank van der Meulen, Moritz Schauer:
    # Sticky PDMP samplers for sparse and local inference problems. 2020. [https://arxiv.org/abs/2103.08478].
    SpikeMixture(basemeasure(μ.m), static(1.0), static(1.0))
end

for func in [:strict_logdensityof, :logdensity_def]
    @eval @inline function $func(μ::SpikeMixture, x)
        # NOTE: We could instead write this as 
        # R1 = typeof(log(one(μ.s))) 
        # R2 = typeof(log(one(μ.w))) 

        # which would rely on constant propagation insteadof type inference.
        # We'll try this for now and come back to the question if we see
        # problems.

        R1 = Core.Compiler.return_type(log, Tuple{typeof(μ.s)})
        R2 = Core.Compiler.return_type(log, Tuple{typeof(μ.w)})
        R3 = Core.Compiler.return_type($func, Tuple{typeof(μ.m),typeof(x)})
        R = promote_type(R1, R2, R3)
        if iszero(x)
            return convert(R, log(μ.s))::R
        else
            return convert(R, log(μ.w) + $func(μ.m, x))::R
        end
    end
end

function gentype(μ::SpikeMixture)
    gentype(μ.m)
end

function Base.rand(rng::AbstractRNG, T::Type, μ::SpikeMixture)
    return (rand(rng, T) < μ.w) * rand(rng, T, μ.m)
end

testvalue(::Type{T}, μ::SpikeMixture) where {T} = zero(T)

insupport(μ::SpikeMixture, x) = dynamic(insupport(μ.m, x)) || iszero(x)
