struct RestrictedMeasure{P,M} <: AbstractMeasure
    predicate::P
    base::M
end

@inline mspace_elsize(μ::RestrictedMeasure) = mspace_elsize(μ.base)
@inline mspace_flatsize(μ::RestrictedMeasure) = mspace_flatsize(μ.base)

@inline logdensity_def(d::RestrictedMeasure, x) = logdensity_def(d.base, x)

basemeasure(μ::RestrictedMeasure) = μ.base

insupport(μ::RestrictedMeasure, x) = _insupport_and(μ.predicate(x), insupport(μ.base, x))

function Pretty.quoteof(d::RestrictedMeasure)
    qf = Pretty.quoteof(d.predicate)
    qbase = Pretty.quoteof(d.base)
    :(RestrictedMeasure($qf, $qbase))
end
