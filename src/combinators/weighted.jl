export WeightedMeasure, AbstractWeightedMeasure

"""
    abstract type AbstractWeightedMeasure <: AbstractMeasure

Supertype of measures that are constant multiples of their base measure,
like [`WeightedMeasure`](@ref) and [`NormalizedMeasure`](@ref). They
sample and transport like their base measure.
"""
abstract type AbstractWeightedMeasure <: AbstractMeasure end

# By default the weight for all measure is 1
_logweight(::AbstractMeasure) = 0
_logweight(μ::AbstractWeightedMeasure) = μ.logweight

@inline logdensity_def(d::AbstractWeightedMeasure, x) = _logweight_for(_logweight(d), x)

# Plain floating-point log-weights adopt the number type of the variate,
# log-weights that carry more information (dual numbers, traced values)
# promote as usual:
@inline _logweight_for(w, x) = w
@inline _logweight_for(w::Union{AbstractFloat,StaticFloat64}, x) = _logd_numtype(x)(dynamic(w))

# The weight-shifted density of a support-safe base density is support-safe,
# no explicit support check required:
@inline function logdensityof_impl(d::AbstractWeightedMeasure, x)
    _add_logweight(_logweight_for(_logweight(d), x), logdensityof_impl(basemeasure(d), x))
end

@inline function batched_logdensityof_impl(d::AbstractWeightedMeasure, X)
    _lazy_add_logweight(_logweight_for(_logweight(d), X), batched_logdensityof_impl(basemeasure(d), X))
end

# A zero weight gives zero densities, also where the base density is
# infinite:
@inline _add_logweight(ℓ_w, ℓ) = (y = ℓ_w + ℓ; ifelse(isneginf(ℓ_w), oftype(y, ℓ_w), y))
@inline _lazy_add_logweight(ℓ_w, ℓ::Number) = _add_logweight(ℓ_w, ℓ)
@inline _lazy_add_logweight(ℓ_w, ℓ) = Broadcast.instantiate(Broadcast.broadcasted(_add_logweight, ℓ_w, ℓ))
@inline function batched_logdensity_def(d::AbstractWeightedMeasure, X)
    _lazy_add(_logweight_for(_logweight(d), X), _zero_logd_batch(X, _static_ndims(basemeasure(d))))
end

@inline rand_impl(ctx::GenContext, μ::AbstractWeightedMeasure) = rand_impl(ctx, basemeasure(μ))
@inline batched_rand_impl(ctx::GenContext, μ::AbstractWeightedMeasure, sz::SizeLike) =
    batched_rand_impl(ctx, basemeasure(μ), sz)

testvalue(::Type{T}, μ::AbstractWeightedMeasure) where {T} = testvalue(T, basemeasure(μ))

###############################################################################

"""
    struct WeightedMeasure{R,M} <: AbstractWeightedMeasure
        logweight :: R
        base :: M
    end
"""
struct WeightedMeasure{R,M} <: AbstractWeightedMeasure
    logweight::R
    base::M
end

@inline mspace_elsize(μ::WeightedMeasure) = mspace_elsize(μ.base)
@inline mspace_flatsize(μ::WeightedMeasure) = mspace_flatsize(μ.base)
@inline mspace_flatsize(::Type{<:WeightedMeasure{<:Any,M}}) where {M} = mspace_flatsize(M)
@inline mspace_ndims(::Type{<:WeightedMeasure{<:Any,M}}) where {M} = mspace_ndims(M)
@inline fixed_stream_size(::Type{<:WeightedMeasure{<:Any,M}}) where {M} = fixed_stream_size(M)

massof(w::WeightedMeasure) = exp(w.logweight) * massof(w.base)

basemeasure(μ::AbstractWeightedMeasure) = μ.base

function Pretty.tile(d::WeightedMeasure)
    weight = round(dynamic(exp(d.logweight)), sigdigits = 4)
    Pretty.pair_layout(Pretty.tile(weight), Pretty.tile(d.base), sep = " * ")
end

function Base.:*(k::T, m::AbstractMeasure) where {T<:Number}
    logk = log(k)
    return weightedmeasure(logk, m)
end

Base.:*(m::AbstractMeasure, k::Number) = k * m

gentype(μ::WeightedMeasure) = gentype(μ.base)

insupportof(μ::WeightedMeasure, x) = insupportof(μ.base, x)

# The Stieltjes measure function scales with the weight:
smf(μ::AbstractWeightedMeasure, x) = _scaled_smf(exp(_logweight(μ)), smf(basemeasure(μ), x))
invsmf(μ::AbstractWeightedMeasure, p) = invsmf(basemeasure(μ), p / exp(_logweight(μ)))

@inline _scaled_smf(w, F) = w * F
@inline _scaled_smf(w, F::NoSMF) = F

# Weighted measures transport like their base:
@inline _unweighted(μ::AbstractWeightedMeasure) = _unweighted(basemeasure(μ))
@inline transport_to_std(::Type{S}, μ::AbstractWeightedMeasure, x) where {S<:StdMeasure} =
    transport_to_std(S, basemeasure(μ), x)
@inline transport_from_std(::Type{S}, μ::AbstractWeightedMeasure, z) where {S<:StdMeasure} =
    transport_from_std(S, basemeasure(μ), z)
@inline transport_to_std_with_rest(::Type{S}, μ::AbstractWeightedMeasure, x::AbstractVector) where {S<:StdMeasure} =
    transport_to_std_with_rest(S, basemeasure(μ), x)
@inline transport_to_std_with_rest(::Type{S}, μ::AbstractWeightedMeasure, x::NamedTuple) where {S<:StdMeasure} =
    transport_to_std_with_rest(S, basemeasure(μ), x)
@inline transport_from_std_with_rest(::Type{S}, μ::AbstractWeightedMeasure, z::AbstractVector) where {S<:StdMeasure} =
    transport_from_std_with_rest(S, basemeasure(μ), z)

@inline batched_transport_to_std(::Type{S}, μ::AbstractWeightedMeasure, X) where {S<:StdMeasure} =
    batched_transport_to_std(S, basemeasure(μ), X)
@inline batched_transport_from_std(::Type{S}, μ::AbstractWeightedMeasure, Z::AbstractArray) where {S<:StdMeasure} =
    batched_transport_from_std(S, basemeasure(μ), Z)
@inline batched_transport_to_std_with_rest(::Type{S}, μ::AbstractWeightedMeasure, X::AbstractArray, sz::SizeLike) where {S<:StdMeasure} =
    batched_transport_to_std_with_rest(S, basemeasure(μ), X, sz)
@inline batched_transport_from_std_with_rest(::Type{S}, μ::AbstractWeightedMeasure, Z::AbstractArray, sz::SizeLike) where {S<:StdMeasure} =
    batched_transport_from_std_with_rest(S, basemeasure(μ), Z, sz)

Adapt.adapt_structure(to, μ::WeightedMeasure) = WeightedMeasure(μ.logweight, Adapt.adapt(to, μ.base))
