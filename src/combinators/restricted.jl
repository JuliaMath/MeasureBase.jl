struct RestrictedMeasure{P,M} <: AbstractMeasure
    predicate::P
    base::M
end

@inline mspace_elsize(μ::RestrictedMeasure) = mspace_elsize(μ.base)
@inline mspace_flatsize(μ::RestrictedMeasure) = mspace_flatsize(μ.base)
@inline mspace_flatsize(::Type{<:RestrictedMeasure{<:Any,M}}) where {M} = mspace_flatsize(M)
@inline mspace_ndims(::Type{<:RestrictedMeasure{<:Any,M}}) where {M} = mspace_ndims(M)
@inline stream_ndims(::Type{<:RestrictedMeasure{<:Any,M}}) where {M} = stream_ndims(M)
@inline fixed_stream_size(::Type{<:RestrictedMeasure{<:Any,M}}) where {M} = fixed_stream_size(M)

# The density relative to the base is one on the restriction, the support
# is handled by `insupportof`:
@inline logdensity_def(::RestrictedMeasure, x) = zero(_logd_numtype(x))

@inline function logdensityof_impl(μ::RestrictedMeasure, x)
    _select(μ.predicate(x), () -> logdensityof(μ.base, x), () -> _neg_inf_logd(x))
end

# Device batches combine a slice-wise predicate with the bulk base density.
# Both must be defined on every variate in the batch.
@inline function batched_logdensityof_impl(μ::RestrictedMeasure, X::AbstractArray{<:Number})
    k = _stream_ndims(μ)
    if k isa StaticInteger && 0 < k < ndims(X) && !_is_host_array(X)
        ℓ = logdensities(μ.base, X)
        if ℓ isa AbstractArray && !_is_host_array(ℓ)
            mask = _broadcast_variate_slices(μ.predicate, X, ℓ)
            return Broadcast.instantiate(Broadcast.broadcasted(_checksupport, mask, ℓ))
        end
    end
    _default_batched_kernel(logdensityof_impl, μ, X, k)
end

basemeasure(μ::RestrictedMeasure) = μ.base
localmeasure(μ::RestrictedMeasure, x) = RestrictedMeasure(μ.predicate, localmeasure(μ.base, x))

insupportof(μ::RestrictedMeasure, x) = _insupport_and(μ.predicate(x), insupportof(μ.base, x))
@inline _relative_support(μ::RestrictedMeasure, x) = μ.predicate(x) & _relative_support(μ.base, x)

Adapt.adapt_structure(to, μ::RestrictedMeasure) = RestrictedMeasure(Adapt.adapt(to, μ.predicate), Adapt.adapt(to, μ.base))

function Pretty.quoteof(d::RestrictedMeasure)
    qf = Pretty.quoteof(d.predicate)
    qbase = Pretty.quoteof(d.base)
    :(RestrictedMeasure($qf, $qbase))
end
