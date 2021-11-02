
###############################################################################
# Affine

affine(f::AffineTransform, μ::AbstractMeasure) = Affine(f, μ)

affine(nt::NamedTuple, μ::AbstractMeasure) = affine(AffineTransform(nt), μ)

affine(f) = μ -> affine(f, μ)

function affine(f::AffineTransform, parent::WeightedMeasure)
    WeightedMeasure(parent.logweight, affine(f, parent.base))
end

function affine(f::AffineTransform, parent::FactoredBase)
    constℓ = parent.constℓ
    varℓ = parent.varℓ
    # Avoid transforming `inbounds`, which is expensive
    base = affine(f, restrict(parent.inbounds, parent.base))
    FactoredBase(Returns(true), constℓ, varℓ, base)
end


###############################################################################
# Half

half(μ::AbstractMeasure) = Half(μ)

###############################################################################
# PointwiseProductMeasure

pointwiseproduct(μ::AbstractMeasure...) = PointwiseProductMeasure(μ)

function pointwiseproduct(
    μ::PointwiseProductMeasure{X},
    ν::PointwiseProductMeasure{Y},
) where {X,Y}
    data = (μ.data..., ν.data...)
    pointwiseproduct(data...)
end

function pointwiseproduct(μ::AbstractMeasure, ν::PointwiseProductMeasure)
    data = (μ, ν.data...)
    return pointwiseproduct(data...)
end

function pointwiseproduct(μ::PointwiseProductMeasure, ν::N) where {N<:AbstractMeasure}
    data = (μ.data..., ν)
    return pointwiseproduct(data...)
end

function pointwiseproduct(μ::AbstractMeasure, ℓ::Likelihood)
    data = (μ, ℓ)
    return PointwiseProductMeasure(data)
end

###############################################################################
# PowerMeaure

function powermeasure(μ::M, dims::NTuple{N,I}) where {M<:AbstractMeasure,N,I}
    productmeasure(Returns(μ), identity, CartesianIndices(dims))
end

function powermeasure(μ::M, dims::Tuple{I}) where {M<:AbstractMeasure,N,I}
    productmeasure(Returns(μ), identity, Base.OneTo(first(dims)))
end

function powermeasure(μ::WeightedMeasure, dims::NTuple{N,I}) where {N,I}
    k = prod(dims) * μ.logweight
    return weightedmeasure(k, μ.base^dims)
end

###############################################################################
# ProductMeasure

productmeasure(f, ops, pars) = ProductMeasure(kernel(f, ops), pars)

productmeasure(f, pars) = productmeasure(f, identity, pars)

productmeasure(μs::Tuple) = TupleProductMeasure(μs)

productmeasure(f::Returns, ops, pars) = ProductMeasure(kernel(f, identity), pars)

productmeasure(k::Kernel, pars) = productmeasure(k.f, k.ops, pars)

productmeasure(nt::NamedTuple) = productmeasure(identity, nt)

function productmeasure(f::Returns{FB}, ops, pars) where {FB<:FactoredBase}
    fb = f.value
    dims = size(pars)
    n = prod(dims)
    inbounds(x) = all(fb.inbounds, x)
    constℓ = n * fb.constℓ
    varℓ() = n * fb.varℓ()
    base = fb.base^dims
    FactoredBase(inbounds, constℓ, varℓ, base)
end

function productmeasure(f::Returns{W}, ::typeof(identity), pars) where {W<:WeightedMeasure}
    ℓ = f.value.logweight
    base = f.value.base
    newbase = productmeasure(Returns(base), identity, pars)
    weightedmeasure(length(pars) * ℓ, newbase)
end



###############################################################################
# RestrictedMeasure

restrict(f, b) = RestrictedMeasure(f, b)

###############################################################################
# SuperpositionMeasure

superpose(a::AbstractArray) = SuperpositionMeasure(a)

superpose(t::Tuple) = SuperpositionMeasure(t)
superpose(nt::NamedTuple) = SuperpositionMeasure(nt)

function superpose(μ::AbstractMeasure, ν::AbstractMeasure)
    components = (μ, ν)
    superpose(components)
end

###############################################################################
# WeightedMeasure

function weightedmeasure(ℓ::R, b::M) where {R,M}
    WeightedMeasure{R,M}(ℓ, b)
end

function weightedmeasure(ℓ, b::WeightedMeasure)
    weightedmeasure(ℓ + b.logweight, b.base)
end

###############################################################################
# Kernel

kernel(μ, ops...) = Kernel(μ, ops)
kernel(μ, op) = Kernel(μ, op)

# kernel(Normal(μ=2))
function kernel(μ::P) where {P<:AbstractMeasure}
    (f, ops) = kernelfactor(μ)
    kernel(f, ops)
end

# kernel(Normal{(:μ,), Tuple{Int64}})
function kernel(::Type{P}) where {P<:AbstractMeasure}
    (f, ops) = kernelfactor(P)
    kernel(f, ops)
end

# kernel(::Type{P}, op::O) where {O, N, P<:ParameterizedMeasure{N}} = kernel{constructorof(P),O}(op)

function kernel(::Type{M}; ops...) where {M}
    nt = NamedTuple(ops)
    kernel(M, nt)
end

kernel(f::Returns, op::typeof(identity)) = Kernel(f, op)
kernel(f::Returns, op) = kernel(f, identity)
