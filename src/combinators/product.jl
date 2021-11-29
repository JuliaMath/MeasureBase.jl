export ProductMeasure

using MappedArrays
using Base: @propagate_inbounds
import Base
using FillArrays

abstract type AbstractProductMeasure <: AbstractMeasure end

struct ProductMeasure{M,K<:AbstractKernel,I} <: AbstractProductMeasure
    f::K
    xs::I

    function ProductMeasure(f::K, xs::I) where {K,I}
        μ₁ = f(first(xs))
        M = typeof(μ₁)
        new{M,K,I}(f, xs)
    end
end

# TODO: Test for equality without traversal, probably by first converting to a
# canonical form
function Base.:(==)(a::ProductMeasure, b::ProductMeasure)
    all(zip(a.xs, b.xs)) do (aᵢ, bᵢ)
        a.f(aᵢ) == b.f(bᵢ)
    end
end

Base.size(μ::ProductMeasure) = size(marginals(μ))

Base.length(μ::ProductMeasure) = length(marginals(μ))

basemeasure(d::ProductMeasure) = productmeasure(basekernel(d.f), d.xs)

basemeasure_depth(μ::ProductMeasure) = basemeasure_depth(first(marginals(μ)))

# TODO: Do we need these methods?
# basemeasure(d::ProductMeasure) = ProductMeasure(basemeasure ∘ d.f, d.xs)
# basemeasure(d::ProductMeasure{typeof(identity)}) = ProductMeasure(identity, map(basemeasure, d.xs))
# basemeasure(d::ProductMeasure{typeof(identity), <:FillArrays.Fill}) = ProductMeasure(identity, map(basemeasure, d.xs))

export marginals

function marginals(d::ProductMeasure{M,K,I}) where {M,K,I}
    _marginals(d, isiterable(I))
end

function _marginals(d::ProductMeasure, ::Iterable)
    return (d.f(i) for i in d.xs)
end

testvalue(d::AbstractProductMeasure) = map(testvalue, marginals(d))

function Base.rand(rng::AbstractRNG, ::Type{T}, d::ProductMeasure) where {T}
    _rand(rng, T, d, marginals(d))
end

function _rand(rng::AbstractRNG, ::Type{T}, d::ProductMeasure, mar) where {T}
    (rand(rng, T, m) for m in mar)
end

###############################################################################
# I <: Tuple

struct TupleProductMeasure{T} <: AbstractProductMeasure
    components::T
end

export ⊗
⊗(μs::AbstractMeasure...) = productmeasure(μs)

marginals(d::TupleProductMeasure{T}) where {F,T<:Tuple} = d.components

@inline function logdensity_def(d::TupleProductMeasure, x::Tuple) where {T<:Tuple}
    mapreduce(logdensity_def, +, d.components, x)
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::TupleProductMeasure) where {T}
    rand.(d.components)
end

###############################################################################
# I <: AbstractArray

marginals(d::ProductMeasure{M,K,A}) where {M,K,A<:AbstractArray} = mappedarray(d.f, d.xs)

@inline function logdensity_def(d::ProductMeasure, x)
    mapreduce(logdensity_def, +, marginals(d), x)
end

@inline function logdensity_def(d::ProductMeasure{<:Returns}, x)
    sum(x -> logdensity_def(d.f.f.value, x), x)
end

# function Base.rand(rng::AbstractRNG, ::Type{T}, d::ProductMeasure{M,K,I}) where {T,F,I<:CartesianIndices}

# end

###############################################################################
# I <: Base.Generator

export rand!
using Random: rand!, GLOBAL_RNG, AbstractRNG


function logdensity_def(μ::AbstractProductMeasure, x)
    mapreduce(+, μ.xs, x) do (j,x)
        logdensity_def(μ.f(j), x)
    end
end

@propagate_inbounds function Random.rand!(
    rng::AbstractRNG,
    d::ProductMeasure,
    x::AbstractArray,
)
    # TODO: Generalize this
    T = Float64
    for (j, m) in zip(eachindex(x), marginals(d))
        @inbounds x[j] = rand(rng, T, m)
    end
    return x
end

export rand!
using Random: rand!, GLOBAL_RNG, AbstractRNG

function _rand(rng::AbstractRNG, ::Type{T}, d::ProductMeasure, mar::AbstractArray) where {T}
    elT = typeof(rand(rng, T, first(mar)))

    sz = size(mar)
    x = Array{elT,length(sz)}(undef, sz)
    rand!(rng, d, x)
end

# TODO: 
# function Base.rand(rng::AbstractRNG, d::ProductMeasure)
#     return rand(rng, gentype(d), d)
# end

# function Base.rand(T::Type, d::ProductMeasure)
#     return rand(Random.GLOBAL_RNG, T, d)
# end

# function Base.rand(d::ProductMeasure)
#     T = gentype(d)
#     return rand(Random.GLOBAL_RNG, T, d)
# end

function gentype(d::ProductMeasure{A}) where {T,N,A<:AbstractArray{T,N}}
    S = @inbounds gentype(marginals(d)[1])
    Array{S,N}
end

function gentype(d::ProductMeasure{<:Tuple})
    Tuple{gentype.(marginals(d))...}
end

# function logdensity_def(μ::ProductMeasure{Aμ}, x::Ax) where {Aμ <: MappedArray, Ax <: AbstractArray}
#     μ.data
# end

# # FIXME
# function ConstructionBase.constructorof(::Type{P}) where {K,I,P<:ProductMeasure{M,K,I}}
#     p -> productmeasure(d.f, p)
# end

# function Accessors.set(d::ProductMeasure{N}, ::typeof(params), p) where {N}
#     setproperties(d, NamedTuple{N}(p...))
# end

# function Accessors.set(d::ProductMeasure{F,T}, ::typeof(params), p::Tuple) where {F, T<:Tuple}
#     set.(marginals(d), params, p)
# end

# function logdensity_def(μ::ProductMeasure, ν::ProductMeasure, x)
#     sum(zip(marginals(μ), marginals(ν), x)) do μ_ν_x
#         logdensity_def(μ_ν_x...)
#     end
# end

function kernelfactor(μ::ProductMeasure{M,K,<:Fill}) where {M,K}
    k = kernel(first(marginals(μ)))
    (p -> k.f(p)^size(μ), k.param_maps)
end

# function kernelfactor(μ::ProductMeasure{K,A}) where {K,A<:AbstractArray}
#     (p -> set.(marginals(μ), params, p), μ.xs)
# end
