
###############################################################################
# Half

half(μ::AbstractMeasure) = Half(μ)

###############################################################################
# PointwiseProductMeasure

function pointwiseproduct(μ::AbstractMeasure, ℓ::Likelihood)
    return PointwiseProductMeasure(μ, ℓ)
end

###############################################################################
# PowerMeaure

function powermeasure(μ::WeightedMeasure, dims::NTuple{N,I}) where {N,I}
    k = mapreduce(length, *, dims) * μ.logweight
    return weightedmeasure(k, μ.base^dims)
end

###############################################################################
# ProductMeasure

productmeasure(mar::Fill) = powermeasure(mar.value, mar.axes)

function productmeasure(mar::ReadonlyMappedArray{T, N, A, Returns{M}}) where {T,N,A,M}
    return powermeasure(mar.f.value, axes(mar.data))
end

productmeasure(mar::Base.Generator) = ProductMeasure(mar)
productmeasure(mar::AbstractArray) = ProductMeasure(mar)


# TODO: Make this static when its length is static
@inline function productmeasure(mar::AbstractArray{WeightedMeasure{StaticFloat64{W}, M}}) where {W,M}
    return weightedmeasure(W * length(mar), productmeasure(map(basemeasure, mar)))
end

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

function superpose(μ::T, ν::T) where {T}
    if μ==ν
        return weightedmeasure(static(logtwo), μ)
    else
        return superpose((μ, ν))
    end
end

function superpose(μ, ν)
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
    kleisli(M)
end

function kleisli(d::PowerMeasure)
    Base.Fix2(powermeasure, d.axes) ∘ kleisli(d.parent)
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