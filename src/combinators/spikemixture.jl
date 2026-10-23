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

for (head, func) in [(:logdensityof_impl, :logdensityof), (:logdensity_def, :logdensity_def)]
    @eval @inline function $head(μ::SpikeMixture, x)
        ℓ_spike = dynamic(log(μ.s))
        ℓ_parent = dynamic(log(μ.w)) + dynamic($func(μ.m, x))
        ifelse(iszero(x), oftype(ℓ_parent, ℓ_spike), ℓ_parent)
    end
end

massof(μ::SpikeMixture) = μ.w * massof(μ.m) + μ.s

function gentype(μ::SpikeMixture)
    gentype(μ.m)
end

# The probability of drawing from the parent rather than the spike:
function _parent_prob(μ::SpikeMixture)
    mass_parent = μ.w * massof(μ.m)
    mass_parent isa AbstractUnknownMass &&
        throw(ArgumentError("Cannot sample from a spike mixture with a parent of unknown mass"))
    return dynamic(mass_parent / (mass_parent + μ.s))
end

function rand_impl(ctx::GenContext, μ::SpikeMixture)
    x = rand_impl(ctx, μ.m)
    return ifelse(rand(ctx) < _parent_prob(μ), x, zero(x))
end

function batched_rand_impl(ctx::GenContext, μ::SpikeMixture, sz::SizeLike)
    _spike_batched_rand(ctx, μ, sz, _stream_ndims(μ.m))
end
function _spike_batched_rand(ctx::GenContext, μ::SpikeMixture, sz::SizeLike, k::StaticInteger)
    X = batched_rand_impl(ctx, μ.m, sz)
    return ifelse.(_batch_mask(_rand_bulk(ctx, sz) .< _parent_prob(μ), k), X, zero(eltype(X)))
end
_spike_batched_rand(ctx::GenContext, μ::SpikeMixture, sz::SizeLike, ::NoMSpaceElementSize) =
    _batched_rand_pointwise(ctx, μ, sz)

testvalue(::Type{T}, μ::SpikeMixture) where {T} = zero(T)

insupportof(μ::SpikeMixture, x) = _insupport_mask(insupportof(μ.m, x)) | iszero(x)


@inline mspace_flatsize(μ::SpikeMixture) = mspace_flatsize(typeof(μ))
@inline mspace_flatsize(::Type{MU}) where {M,MU<:SpikeMixture{M}} = _scalar_or_unknown(MU, mspace_flatsize(M))
