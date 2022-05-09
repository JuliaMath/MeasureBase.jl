
###############################################################################
# Half

half(μ::AbstractMeasure) = Half(μ)

###############################################################################
# PointwiseProductMeasure

function pointwiseproduct(μ::AbstractMeasure, ℓ::Likelihood)
    T = Core.Compiler.return_type(ℓ.k, Tuple{gentype(μ)})
    return pointwiseproduct(T, μ, ℓ)
end

function pointwiseproduct(::Type{T}, μ::AbstractMeasure, ℓ::Likelihood) where {T}
    return PointwiseProductMeasure(μ, ℓ)
end

###############################################################################
# PowerMeaure

function powermeasure(μ::WeightedMeasure, dims::NTuple{N,I}) where {N,I<:AbstractArray}
    k = mapreduce(length, *, dims) * μ.logweight
    return weightedmeasure(k, μ.base^dims)
end

function powermeasure(μ::WeightedMeasure, dims::NTuple{N,I}) where {N,I}
    k = prod(dims) * μ.logweight
    return weightedmeasure(k, μ.base^dims)
end

###############################################################################
# ProductMeasure

productmeasure(mar::Fill) = powermeasure(mar.value, mar.axes)

function productmeasure(mar::ReadonlyMappedArray{T,N,A,Returns{M}}) where {T,N,A,M}
    return powermeasure(mar.f.value, axes(mar.data))
end

productmeasure(mar::Base.Generator) = ProductMeasure(mar)
productmeasure(mar::AbstractArray) = ProductMeasure(mar)

# TODO: Make this static when its length is static
@inline function productmeasure(
    mar::AbstractArray{WeightedMeasure{StaticFloat64{W},M}},
) where {W,M}
    return weightedmeasure(W * length(mar), productmeasure(map(basemeasure, mar)))
end

productmeasure(nt::NamedTuple) = ProductMeasure(nt)
productmeasure(tup::Tuple) = ProductMeasure(tup)

productmeasure(f, param_maps, pars) = ProductMeasure(kernel(f, param_maps), pars)

function productmeasure(k::ParameterizedTransitionKernel, pars)
    productmeasure(k.f, k.param_maps, pars)
end

function productmeasure(f::Returns{W}, ::typeof(identity), pars) where {W<:WeightedMeasure}
    ℓ = _logweight(f.value)
    base = basemeasure(f.value)
    newbase = productmeasure(Returns(base), identity, pars)
    weightedmeasure(length(pars) * ℓ, newbase)
end

###############################################################################
# RestrictedMeasure
export restrict

restrict(f, b) = RestrictedMeasure(f, b)

###############################################################################
# SuperpositionMeasure

superpose(a::AbstractArray) = SuperpositionMeasure(a)

superpose(t::Tuple) = SuperpositionMeasure(t)
superpose(nt::NamedTuple) = SuperpositionMeasure(nt)

function superpose(μ::T, ν::T) where {T<:AbstractMeasure}
    if μ == ν
        return weightedmeasure(static(logtwo), μ)
    else
        return superpose((μ, ν))
    end
end

function superpose(μ::AbstractMeasure, μs...)
    if all(==(μ), μs)
        return weightedmeasure(log(length(μs) + 1), μ)
    else
        return superpose((μ, μs...))
    end
end

add_measures(μs::AbstractVector, νs) = push!(μs, νs...)
add_measures(μs::Tuple, νs) = (μs..., νs...)

function superpose(μ::SuperpositionMeasure, μs...)
    SuperpositionMeasure(add_measures(μ.components, μs))
end

superpose(μ::SuperpositionMeasure) = μ

###############################################################################
# WeightedMeasure

function weightedmeasure(ℓ::R, b::M) where {R,M}
    WeightedMeasure{R,M}(ℓ, b)
end

function weightedmeasure(ℓ, b::WeightedMeasure)
    weightedmeasure(ℓ + _logweight(b), b.base)
end

###############################################################################
# TransitionKernel

kernel(f, pars::NamedTuple) = ParameterizedTransitionKernel(f, pars)

# kernel(Normal(μ=2))
function kernel(μ::M) where {M<:ParameterizedMeasure}
    kernel(M)
end

function kernel(d::PowerMeasure)
    Base.Fix2(powermeasure, d.axes) ∘ kernel(d.parent)
end

# kernel(Normal{(:μ,), Tuple{Int64}})
function kernel(::Type{M}) where {M<:AbstractMeasure}
    constructorof(M)
end

# kernel(::Type{P}, op::O) where {O, N, P<:ParameterizedMeasure{N}} = kernel{constructorof(P),O}(op)

function kernel(::Type{M}; param_maps...) where {M}
    nt = NamedTuple(param_maps)
    kernel(M, nt)
end

kernel(k::ParameterizedTransitionKernel) = k
