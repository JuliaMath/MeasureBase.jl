export FactoredBase

struct FactoredBase{R,C,V,B} <: AbstractMeasure
    inbounds::R
    constℓ::C
    varℓ::V
    base::B
end

@inline function logdensity_def(d::FactoredBase, x)
    d.inbounds(x) || return -Inf
    d.constℓ + d.varℓ()
end

function Pretty.tile(fb::FactoredBase)
    result = Pretty.literal("FactoredBase")
    result *= Pretty.list_layout(Pretty.tile.([fb.inbounds, fb.constℓ, fb.varℓ, fb.base]))
    result
end

basemeasure(d::FactoredBase) = d.base
