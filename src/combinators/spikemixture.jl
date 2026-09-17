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

for func in [:logdensityof, :logdensity_def]
    @eval @inline function $func(μ::SpikeMixture, x)
        ℓ_spike = dynamic(log(μ.s))
        ℓ_parent = dynamic(log(μ.w)) + dynamic($func(μ.m, x))
        ifelse(iszero(x), oftype(ℓ_parent, ℓ_spike), ℓ_parent)
    end
end

function gentype(μ::SpikeMixture)
    gentype(μ.m)
end

function rand_impl(ctx::GenContext, μ::SpikeMixture)
    return (rand(get_rng(ctx), get_precision(ctx)) < μ.w) * rand_impl(ctx, μ.m)
end

function batched_rand_impl(ctx::GenContext, μ::SpikeMixture, sz::Dims)
    X = batched_rand_impl(ctx, μ.m, sz)
    return ifelse.(_rand_bulk(ctx, sz) .< μ.w, X, zero(eltype(X)))
end

testvalue(::Type{T}, μ::SpikeMixture) where {T} = zero(T)

insupport(μ::SpikeMixture, x) = _insupport_mask(insupport(μ.m, x)) | iszero(x)


@inline mspace_flatsize(μ::SpikeMixture) = _scalar_or_unknown(mspace_flatsize(μ.m))
@inline mspace_flatsize(::Type{<:SpikeMixture{M}}) where {M} = _scalar_or_unknown(mspace_flatsize(M))
