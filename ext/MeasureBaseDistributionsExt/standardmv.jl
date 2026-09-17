# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).


MeasureBase.getdof(d::AbstractMvNormal) = length(d)
MeasureBase.getdof(m::AsMeasure{<:AbstractMvNormal}) = getdof(m.obj)

@inline MeasureBase.preferred_stdmeasure(::Type{<:AbstractMvNormal}) = StdNormal

_cholesky_L(A) = cholesky(A).L
_cholesky_L(A::Diagonal{<:Real}) = Diagonal(sqrt.(diag(A)))
_cholesky_L(A::PDMats.PDiagMat{<:Real}) = Diagonal(sqrt.(A.diag))
_cholesky_L(A::PDMats.ScalMat{<:Real}) = Diagonal(Fill(sqrt(A.value), A.dim))

function MeasureBase.transport_to_std(::Type{StdNormal}, d::MvNormal, x)
    _cholesky_L(d.Σ) \ (x - d.μ)
end

function MeasureBase.transport_from_std(::Type{StdNormal}, d::MvNormal, z)
    muladd(_cholesky_L(d.Σ), z, d.μ)
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
