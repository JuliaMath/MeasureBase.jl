# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).


MeasureBase.getdof(d::AbstractMvNormal) = length(d)
MeasureBase.getdof(m::AsMeasure{<:AbstractMvNormal}) = getdof(m.obj)

@inline MeasureBase.preferred_stdmeasure(::Type{<:AbstractMvNormal}) = StdNormal

_cholesky_L(A) = _lower_factor(cholesky(A))
# `Cholesky.L` copies the transposed factor and scalar-indexes on GPUs.
_lower_factor(C::Cholesky) = C.uplo === 'L' ? LowerTriangular(C.factors) : UpperTriangular(C.factors)'
_cholesky_L(A::Diagonal{<:Real}) = Diagonal(sqrt.(diag(A)))
_cholesky_L(A::PDMats.PDiagMat{<:Real}) = Diagonal(sqrt.(A.diag))
_cholesky_L(A::PDMats.ScalMat{<:Real}) = Diagonal(Fill(sqrt(A.value), A.dim))

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
