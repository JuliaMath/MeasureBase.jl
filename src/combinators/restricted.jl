struct RestrictedMeasure{P,M} <: AbstractMeasure
    predicate::P
    base::M
end

@inline mspace_elsize(μ::RestrictedMeasure) = mspace_elsize(μ.base)
@inline mspace_flatsize(μ::RestrictedMeasure) = mspace_flatsize(μ.base)
@inline mspace_flatsize(::Type{<:RestrictedMeasure{<:Any,M}}) where {M} = mspace_flatsize(M)
@inline mspace_ndims(::Type{<:RestrictedMeasure{<:Any,M}}) where {M} = mspace_ndims(M)
@inline fixed_stream_size(::Type{<:RestrictedMeasure{<:Any,M}}) where {M} = fixed_stream_size(M)

@inline logdensity_def(d::RestrictedMeasure, x) = logdensity_def(d.base, x)

basemeasure(μ::RestrictedMeasure) = μ.base

insupport(μ::RestrictedMeasure, x) = _insupport_and(μ.predicate(x), insupport(μ.base, x))

function Pretty.quoteof(d::RestrictedMeasure)
    qf = Pretty.quoteof(d.predicate)
    qbase = Pretty.quoteof(d.base)
    :(RestrictedMeasure($qf, $qbase))
end
