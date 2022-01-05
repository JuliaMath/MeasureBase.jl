export ProductMeasure

using MappedArrays
using MappedArrays: ReadonlyMultiMappedArray
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

function Pretty.tile(d::ProductMeasure{T}) where {T<:Tuple}
    Pretty.list_layout(Pretty.tile.([marginals(d)...]), sep=" ⊗ ")
end

# For tuples, `mapreduce` has trouble with type inference
@inline function logdensity_def(d::ProductMeasure{T}, x) where {T<:Tuple}
    ℓs = map(logdensity_def, marginals(d),x)
    sum(ℓs)
end

function basemeasure(μ::ProductMeasure{Base.Generator{I,F}}) where {I,F}
    T = Core.Compiler.return_type(mar.f, Tuple{_eltype(mar.iter)})
    B = Core.Compiler.return_type(basemeasure, Tuple{T})
    _basemeasure(μ, B, static(Base.issingulartype(B)))
end



function basemeasure(μ::ProductMeasure{A}) where {T,A<:AbstractMappedArray{T}}
    B = Core.Compiler.return_type(basemeasure, Tuple{T})
    _basemeasure(μ, B, static(Base.issingulartype(B)))
end

function _basemeasure(μ::ProductMeasure, ::Type{B}, ::True) where {T,B}
    return instance(B) ^ axes(marginals(μ))
end

function _basemeasure(μ::ProductMeasure{A}, ::Type{B}, ::False) where {T,A<:AbstractMappedArray{T},B}
    productmeasure(mappedarray(basemeasure, mar))
end

function _basemeasure(μ::ProductMeasure{Base.Generator{I,F}}, ::Type{B}, ::False) where {I,F,B}
    mar = marginals(μ)
    productmeasure(Base.Generator(basekleisli(mar.f), mar.iter))
end

marginals(μ::ProductMeasure) = μ.marginals


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
