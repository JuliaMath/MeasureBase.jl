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

basemeasure(μ::RestrictedMeasure) = μ.base

insupportof(μ::RestrictedMeasure, x) = _insupport_and(μ.predicate(x), insupportof(μ.base, x))

function Pretty.quoteof(d::RestrictedMeasure)
    qf = Pretty.quoteof(d.predicate)
    qbase = Pretty.quoteof(d.base)
    :(RestrictedMeasure($qf, $qbase))
end
