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
    mar = marginals(μ)
    T = Core.Compiler.return_type(mar.f, Tuple{_eltype(mar.iter)})
    if static_hasmethod(tbasemeasure_type, Tuple{Type{T}})
        B = tbasemeasure_type(T)
        if Base.issingletontype(B)
            b = instance(B)::B
            return b ^ axes(mar)
        end
    end

    return productmeasure(Base.Generator(basekleisli(mar.f), mar.iter))
end

function basemeasure(μ::ProductMeasure{A}) where {T,A<:AbstractMappedArray{T}}
    mar = marginals(μ)
    if static_hasmethod(tbasemeasure_type, Tuple{Type{T}})
        B = tbasemeasure_type(T)
        if Base.issingletontype(B)
            b = instance(B)
            mar = marginals(μ)
            return b ^ axes(marginals(μ))
        end
    end
    return productmeasure(mappedarray(basemeasure, mar))
end

marginals(μ::ProductMeasure) = μ.marginals

@inline function tbasemeasure_type(::Type{ProductMeasure{M}}) where {T,N,A,F,M<:ReadonlyMultiMappedArray{T,N,A,F}}
    if static_hasmethod(tbasemeasure_type, Tuple{Type{T}})
        B = tbasemeasure_type(T)
        if Base.issingletontype(B)
            return PowerMeasure{B, A}
        end
    end
    tmap(tbasemeasure_type, ProductMeasure{M})
end

@inline function tbasemeasure_type(::Type{P}) where {T,N,A,F,M<:ReadonlyMappedArray{T,N,A,F},P<:ProductMeasure{M}}
    if static_hasmethod(tbasemeasure_type, Tuple{Type{T}})
        B = tbasemeasure_type(T)
        if Base.issingletontype(B)
            return PowerMeasure{B, Tuple{A}}
        else
            tmap(tbasemeasure_type, P)
        end
    else
        @error "tbasemeasure_type(::Type{$P})"
    end
end

@inline function tbasemeasure_type(::Type{P}) where {P<:ProductMeasure}
    tmap(tbasemeasure_type, P)
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
