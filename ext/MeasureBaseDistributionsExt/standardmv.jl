# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).


MeasureBase.getdof(d::AbstractMvNormal) = length(d)
MeasureBase.getdof(m::AsMeasure{<:AbstractMvNormal}) = getdof(m.obj)

MeasureBase.transport_origin(ν::MvNormal) = StandardDist{Normal}(length(ν))

_cholesky_L(A) = cholesky(A).L
_cholesky_L(A::Diagonal{<:Real}) = Diagonal(sqrt.(diag(A)))
_cholesky_L(A::PDMats.PDiagMat{<:Real}) = Diagonal(sqrt.(A.diag))
_cholesky_L(A::PDMats.ScalMat{<:Real}) = Diagonal(Fill(sqrt(A.value), A.dim))

function MeasureBase.from_origin(ν::MvNormal, x)
    A = _cholesky_L(ν.Σ)
    b = ν.μ
    muladd(A, x, b)
end

function MeasureBase.to_origin(ν::MvNormal, y)
    A = _cholesky_L(ν.Σ)
    b = ν.μ
    A \ (y - b)
end


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
