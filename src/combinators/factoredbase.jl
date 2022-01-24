export FactoredBase

struct FactoredBase{C,V,B} <: AbstractMeasure
    constℓ::C
    varℓ::V
    base::B
end

@inline function logdensity_def(d::FactoredBase, x)
    d.constℓ + d.varℓ()
end

function Pretty.tile(fb::FactoredBase)
    result = Pretty.literal("FactoredBase")
    result *= Pretty.list_layout(Pretty.tile.([fb.constℓ, fb.varℓ, fb.base]))
    result
end

basemeasure(d::FactoredBase) = d.base
