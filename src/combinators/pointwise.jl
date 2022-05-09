export ⊙

struct PointwiseProductMeasure{P,L} <: AbstractMeasure
    prior::P
    likelihood::L
end



iterate(p::PointwiseProductMeasure, i=1) = iterate((p.prior, p.likelihood), i)

function Pretty.tile(d::PointwiseProductMeasure)
    Pretty.pair_layout(Pretty.tile(d.prior), Pretty.tile(d.likelihood), sep=" ⊙ ")
end

⊙(μ, ℓ) = pointwiseproduct(μ, ℓ)

@inline function logdensity_def(d::PointwiseProductMeasure, p)
    μ, ℓ = d
    logdensityof(ℓ.k(p), ℓ.x)
end

function gentype(d::PointwiseProductMeasure)
    gentype(d.prior)
end

@inbounds function insupport(d::PointwiseProductMeasure, p) 
    μ, ℓ = d
    insupport(μ, p) && insupport(ℓ.k(p), ℓ.x)
end

basemeasure(d::PointwiseProductMeasure, x) = d.prior

basemeasure(d::PointwiseProductMeasure) = basemeasure(d.prior)
