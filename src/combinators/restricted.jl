struct RestrictedMeasure{F,M} <: AbstractMeasure
    f::F
    base::M
end

@inline function logdensity(d::RestrictedMeasure, x)
    d.f(x) || return -Inf
    return 0.0
end

function density(d::RestrictedMeasure, x)
    d.f(x) || return 0.0
    return 1.0
end

basemeasure(μ::RestrictedMeasure) = μ.base
