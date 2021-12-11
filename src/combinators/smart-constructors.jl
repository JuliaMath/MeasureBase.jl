
###############################################################################
# Half

half(μ::AbstractMeasure) = Half(μ)

###############################################################################
# PointwiseProductMeasure

pointwiseproduct(μ::AbstractMeasure...) = PointwiseProductMeasure(μ)

function pointwiseproduct(μ::AbstractMeasure, ℓ::Likelihood)
    data = (μ, ℓ)
    return PointwiseProductMeasure(data)
end

###############################################################################
# PowerMeaure

function powermeasure(μ::WeightedMeasure, dims::NTuple{N,I}) where {N,I}
    k = mapreduce(length, *, dims) * μ.logweight
    return weightedmeasure(k, μ.base^dims)
end

###############################################################################
# ProductMeasure

productmeasure(data::AbstractArray) = ProductMeasure(data)
productmeasure(nt::NamedTuple) = ProductMeasure(nt)
productmeasure(tup::Tuple) = ProductMeasure(tup)

productmeasure(f, param_maps, pars) = ProductMeasure(kleisli(f, param_maps), pars)

productmeasure(k::ParameterizedKleisli, pars) = productmeasure(k.f, k.param_maps, pars)


function productmeasure(f::Returns{FB}, param_maps, pars) where {FB<:FactoredBase}
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
# Kleisli

kleisli(μ, op1, op2, param_maps...) = ParameterizedKleisli(μ, op1, op2, param_maps...)

# kleisli(Normal(μ=2))
function kleisli(μ::M) where {M<:AbstractMeasure}
    constructorof(M)
end

# kleisli(Normal{(:μ,), Tuple{Int64}})
function kleisli(::Type{M}) where {M<:AbstractMeasure}
    constructorof(M)
end

# kleisli(::Type{P}, op::O) where {O, N, P<:ParameterizedMeasure{N}} = kleisli{constructorof(P),O}(op)

function kleisli(::Type{M}; param_maps...) where {M}
    nt = NamedTuple(param_maps)
    kleisli(M, nt)
end


kleisli(k::ParameterizedKleisli) = k