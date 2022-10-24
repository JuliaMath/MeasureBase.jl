export ⊙

struct PointwiseProductMeasure{P,L} <: AbstractMeasure
    prior::P
    likelihood::L
end

iterate(p::PointwiseProductMeasure, i = 1) = iterate((p.prior, p.likelihood), i)

function Pretty.tile(d::PointwiseProductMeasure)
    Pretty.pair_layout(Pretty.tile(d.prior), Pretty.tile(d.likelihood), sep = " ⊙ ")
end

⊙(prior, ℓ) = pointwiseproduct(prior, ℓ)

@inbounds function insupport(d::PointwiseProductMeasure, p)
    prior, ℓ = d
    istrue(insupport(prior, p)) && istrue(insupport(ℓ, p))
end

@inline function logdensity_def(d::PointwiseProductMeasure, p)
    prior, ℓ = d
    unsafe_logdensityof(ℓ, p)
end

basemeasure(d::PointwiseProductMeasure) = d.prior

function gentype(d::PointwiseProductMeasure)
    gentype(d.prior)
end
