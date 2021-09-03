
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
