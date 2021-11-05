"""
    FactoredBase{R,C,V,B} <: AbstractMeasure

Wrapper around a base measure that stores several operations we need.

# Fields
- `inbounds::R`: function of `x` checking whether it is in the domain of `base`
- `constℓ::C`: constant term in the logdensity of `base` (independent of `x` and `params(base)`)
- `varℓ::V`: variable term in the logdensity of `base` (independent of `x` but not of `params(base)`)
- `base::B`: actual base measure
"""
struct FactoredBase{R,C,V,B} <: AbstractMeasure
    inbounds::R
    constℓ::C
    varℓ::V
    base::B
end

@inline function logdensity(d::FactoredBase, x)
    d.inbounds(x) || return -Inf
    d.constℓ + d.varℓ()
end

basemeasure(d::FactoredBase) = d.base
