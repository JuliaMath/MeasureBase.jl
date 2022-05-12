struct RestrictedMeasure{P,M} <: AbstractMeasure
    predicate::P
    base::M
end

@inline logdensity_def(d::RestrictedMeasure, x) = logdensity_def(d.base, x)

basemeasure(μ::RestrictedMeasure) = μ.base

insupport(μ::RestrictedMeasure, x) = μ.predicate(x) && insupport(μ.base, x)

function Pretty.quoteof(d::RestrictedMeasure)
    qf = Pretty.quoteof(d.predicate)
    qbase = Pretty.quoteof(d.base)
    :(RestrictedMeasure($qf, $qbase))
end
