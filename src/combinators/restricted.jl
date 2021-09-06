struct RestrictedMeasure{F,M} <: AbstractMeasure
    f::F
    base::M
end

function logdensity(d::Restricted, x)
    d.f(x) || return -Inf
end

function density(d::Restricted, x)
    d.f(x) || return 0.0
end

basemeasure(μ::RestrictedMeasure) = μ.base
