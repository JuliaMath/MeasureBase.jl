struct RestrictedMeasure{F,M} <: AbstractMeasure
    f::F
    base::M
end

restrict(f, b) = RestrictedMeasure(f, b)

function logdensity(d::Restricted, x)
    d.f(x) || return -Inf
    logdensity(d.base, x)
end

function density(d::Restricted, x)
    d.f(x) || return 0.0
    density(d.base, x)
end
