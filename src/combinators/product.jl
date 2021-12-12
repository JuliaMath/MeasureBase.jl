export ProductMeasure

using MappedArrays
using Base: @propagate_inbounds
import Base
using FillArrays

abstract type AbstractProductMeasure <: AbstractMeasure end

function Pretty.tile(μ::AbstractProductMeasure)
    result = Pretty.literal("ProductMeasure(")
    result *= Pretty.tile(marginals(μ))
    result *= Pretty.literal(")")
end

export marginals

function Base.:(==)(a::AbstractProductMeasure, b::AbstractProductMeasure)
    marginals(a) == marginals(b)
end
Base.length(μ::AbstractProductMeasure) = length(marginals(μ))
Base.size(μ::AbstractProductMeasure) = size(marginals(μ))
basemeasure(d::AbstractProductMeasure) = productmeasure(map(basemeasure, marginals(d)))

function Base.rand(rng::AbstractRNG, ::Type{T}, d::AbstractProductMeasure) where {T}
    map(marginals(d)) do dⱼ
        rand(rng, T, dⱼ)
    end
end

@inline function logdensity_def(d::AbstractProductMeasure, x)
    mapreduce(logdensity_def, +, marginals(d), x)
end

struct ProductMeasure{M} <: AbstractProductMeasure
    marginals::M
end

@inline function logdensityof(μ::ProductMeasure, x)
    mapreduce(logdensityof, +, marginals(μ), x)
end

marginals(μ::ProductMeasure) = μ.marginals

@generated function tbasemeasure_type(::Type{ProductMeasure{T}}) where {T<:Tuple}
    ProductMeasure{Tuple{map(tbasemeasure_type, T.types)...}}
end

function tbasemeasure_type(::Type{ProductMeasure{A}}) where {M,A<:AbstractArray{M}}
    C = constructorof(A)
    p = Tuple(A.parameters)
    ProductMeasure{C{(tbasemeasure_type(first(p)), Base.tail(p)...)...}}
end

# basemeasure_depth(μ::ProductMeasure) = basemeasure_depth(first(marginals(μ)))

@inline function tbasemeasure_depth(::Type{ProductMeasure{T}}) where {T<:Tuple}
    mapreduce(tbasemeasure_depth, max, Tuple(T.types))
end


testvalue(d::AbstractProductMeasure) = map(testvalue, marginals(d))



export ⊗

"""
    ⊗(μs::AbstractMeasure...)

`⊗` is a binary operator for building product measures. This satisfies the law

```
    basemeasure(μ ⊗ ν) == basemeasure(μ) ⊗ basemeasure(ν)
```
"""
⊗(μs::AbstractMeasure...) = productmeasure(μs)

###############################################################################
# I <: Base.Generator

export rand!
using Random: rand!, GLOBAL_RNG, AbstractRNG

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
