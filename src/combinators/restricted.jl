struct RestrictedMeasure{F,M} <: AbstractMeasure
    f::F
    base::M
end

@inline function logdensity_def(d::RestrictedMeasure, x)
    d.f(x) || return -Inf
    return 0.0
end

function density_def(d::RestrictedMeasure, x)
    d.f(x) || return 0.0
    return 1.0
end

basemeasure(μ::RestrictedMeasure) = μ.base

basemeasure_depth(::Type{RestrictedMeasure{F,M}}) where {F,M} =
    static(1) + basemeasure_depth(M)
