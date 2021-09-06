
###############################################################################
# Affine

affine(nt::NamedTuple, μ::AbstractMeasure) = affine(AffineTransform(nt), μ)

affine(nt::NamedTuple) = μ -> affine(nt, μ)

function affine(f::AffineTransform, parent::WeightedMeasure)
    WeightedMeasure(parent.logweight, affine(f, parent.base))
end

function affine(f::AffineTransform, parent::FactoredBase)
    finv = inv(f)
    inbounds(x) = parent.inbounds(finv(x))
    constℓ = parent.constℓ
    varℓ = parent.varℓ
    base = AffineTransform(f, parent.base)
    FactoredBase(inbounds, constℓ, varℓ, base)
end


###############################################################################
# Half

function half end

###############################################################################
# PointwiseProductMeasure

function pointwiseproduct end

###############################################################################
# PowerMeaure

function powermeasure(μ::M, dims::NTuple{N,I}) where {M<:AbstractMeasure,N,I<:Integer}
    productmeasure(identity, Fill(μ, dims))
end

function powermeasure(μ::WeightedMeasure, dims::NTuple{N,I}) where {N,I<:Integer}
    k = prod(dims) * μ.logweight
    return WeightedMeasure(k, μ.base^dims)
end

###############################################################################
# ProductMeasure

function productmeasure end

productmeasure(nt::NamedTuple) = productmeasure(identity, nt)

###############################################################################
# RestrictedMeasure

restrict(f, b) = RestrictedMeasure(f, b)

###############################################################################
# SuperpositionMeasure

function superpose(μ...)
end


###############################################################################
# WeightedMeasure

function weightedmeasure(ℓ::R, b::M) where {R,M}
    WeightedMeasure{R,M}(ℓ, b)
end 

function weightedmeasure(ℓ, b::WeightedMeasure)
    weightedmeasure(ℓ + b.logweight, b.base)
end 

