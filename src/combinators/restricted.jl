struct RestrictedMeasure{F,M} <: AbstractMeasure
    f::F
    base::M
end

function logdensity(d::RestrictedMeasure, x)
    d.f(x) || return -Inf
end

function density(d::RestrictedMeasure, x)
    d.f(x) || return 0.0
end

basemeasure(μ::RestrictedMeasure) = μ.base
