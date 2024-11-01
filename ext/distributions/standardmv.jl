# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).


MeasureBase.getdof(m::AsMeasure{<:AbstractMvNormal}) = length(m.obj)

MeasureBase.transport_origin(ν::MvNormal) = StandardDist{Normal}(length(ν))

function MeasureBase.from_origin(ν::MvNormal, x)
    A = cholesky(ν.Σ).L
    b = ν.μ
    muladd(A, x, b)
end

function MeasureBase.to_origin(ν::MvNormal, y)
    A = cholesky(ν.Σ).L
    b = ν.μ
    A \ (y - b)
end


AbstractMvNormal
AbstractMvLogNormal

#DirichletMultinomial
#Distributions.AbstractMvLogNormal
#Distributions.AbstractMvTDist
#Distributions.ProductDistribution{1}
#Distributions.ReshapedDistribution{1, S, D} where {S<:ValueSupport, D<:(Distribution{<:ArrayLikeVariate, S})}
#JointOrderStatistics
#Multinomial
#MultivariateMixture (alias for AbstractMixtureModel{ArrayLikeVariate{1}})
#MvLogitNormal
#VonMisesFisher
