export FactoredBase

struct FactoredBase{R,C,V,B} <: AbstractMeasure
    inbounds :: R
    constℓ :: C
    varℓ :: V
    base :: B
end

function logdensity(d::FactoredBase, x)
    d.inbounds(x) || return -Inf
    d.constℓ + d.varℓ
end

basemeasure(d::FactoredBase) = d.base
