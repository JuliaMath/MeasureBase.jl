export WeightedMeasure, AbstractWeightedMeasure

"""
    struct WeightedMeasure{R,M} <: AbstractMeasure
        logweight :: R
        base :: M
    end
"""

abstract type AbstractWeightedMeasure <: AbstractMeasure end

# By default the weight for all measure is 1
_logweight(::AbstractMeasure) = 0

@inline logdensity_def(d::AbstractWeightedMeasure, x) = _logweight_for(d.logweight, x)

# Plain floating-point log-weights adopt the number type of the variate,
# log-weights that carry more information (dual numbers, traced values)
# promote as usual:
@inline _logweight_for(w, x) = w
@inline _logweight_for(w::Union{AbstractFloat,StaticFloat64}, x) = _logd_numtype(x)(dynamic(w))

# The weight-shifted density of a support-safe base density is support-safe,
# no explicit support check required:
@inline function logdensityof_impl(d::AbstractWeightedMeasure, x)
    _logweight_for(d.logweight, x) + logdensityof_impl(basemeasure(d), x)
end

@inline function batched_logdensityof_impl(d::AbstractWeightedMeasure, X)
    _lazy_add(_logweight_for(d.logweight, X), batched_logdensityof_impl(basemeasure(d), X))
end
@inline function batched_logdensity_def(d::AbstractWeightedMeasure, X)
    _lazy_add(_logweight_for(d.logweight, X), _zero_logd_batch(X, mspace_ndims(basemeasure(d))))
end

@inline rand_impl(ctx::GenContext, μ::AbstractWeightedMeasure) = rand_impl(ctx, basemeasure(μ))
@inline batched_rand_impl(ctx::GenContext, μ::AbstractWeightedMeasure, sz::Dims) =
    batched_rand_impl(ctx, basemeasure(μ), sz)

testvalue(::Type{T}, μ::AbstractWeightedMeasure) where {T} = testvalue(T, basemeasure(μ))

###############################################################################

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

_logweight(μ::WeightedMeasure) = μ.logweight
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

insupport(μ::WeightedMeasure, x) = insupport(μ.base, x)

# Weighted measures transport like their base:
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

@inline batched_transport_to_std(::Type{S}, μ::AbstractWeightedMeasure, X::AbstractArray) where {S<:StdMeasure} =
    batched_transport_to_std(S, basemeasure(μ), X)
@inline batched_transport_from_std(::Type{S}, μ::AbstractWeightedMeasure, Z::AbstractArray) where {S<:StdMeasure} =
    batched_transport_from_std(S, basemeasure(μ), Z)
@inline batched_transport_to_std_with_rest(::Type{S}, μ::AbstractWeightedMeasure, X::AbstractArray) where {S<:StdMeasure} =
    batched_transport_to_std_with_rest(S, basemeasure(μ), X)
@inline batched_transport_from_std_with_rest(::Type{S}, μ::AbstractWeightedMeasure, Z::AbstractArray) where {S<:StdMeasure} =
    batched_transport_from_std_with_rest(S, basemeasure(μ), Z)

Adapt.adapt_structure(to, μ::WeightedMeasure) = WeightedMeasure(μ.logweight, Adapt.adapt(to, μ.base))
