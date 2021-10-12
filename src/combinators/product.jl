export ProductMeasure

using MappedArrays
using Base: @propagate_inbounds
import Base
using FillArrays

abstract type AbstractProductMeasure <: AbstractMeasure end


struct ProductMeasure{F,S,I} <: AbstractProductMeasure
    f::Kernel{F,S}
    pars::I
end


# TODO: Test for equality without traversal, probably by first converting to a
# canonical form
function Base.:(==)(a::ProductMeasure, b::ProductMeasure)
    all(zip(a.pars, b.pars)) do (aᵢ, bᵢ)
        a.f(aᵢ) == b.f(bᵢ)
    end
end

Base.size(μ::ProductMeasure) = size(marginals(μ))

Base.length(m::ProductMeasure{T}) where {T} = length(marginals(μ))

basemeasure(d::ProductMeasure) = productmeasure(basekernel(d.f), d.pars)

# TODO: Do we need these methods?
# basemeasure(d::ProductMeasure) = ProductMeasure(basemeasure ∘ d.f, d.pars)
# basemeasure(d::ProductMeasure{typeof(identity)}) = ProductMeasure(identity, map(basemeasure, d.pars))
# basemeasure(d::ProductMeasure{typeof(identity), <:FillArrays.Fill}) = ProductMeasure(identity, map(basemeasure, d.pars))

export marginals

function marginals(d::ProductMeasure{F,S,I}) where {F,S,I}
    _marginals(d, isiterable(I))
end

function _marginals(d::ProductMeasure, ::Iterable)
    return (d.f(i) for i in d.pars)
end

function _marginals(d::ProductMeasure{F,S,I}, ::NonIterable) where {F,S,I}
    error("Type $I is not iterable. Add an `iterate` or `marginals` method to fix.")
end

testvalue(d::ProductMeasure) = map(testvalue, marginals(d))

function Base.show(io::IO, μ::ProductMeasure{F,S,NamedTuple{N,T}}) where {F,S,N,T}
    io = IOContext(io, :compact => true)
    print(io, "Product(",μ.data, ")")
end

function Base.show_unquoted(io::IO, μ::ProductMeasure, indent::Int, prec::Int)
    io = IOContext(io, :compact => true)
    if Base.operator_precedence(:*) ≤ prec
        print(io, "(")
        show(io, μ)
        print(io, ")")
    else
        show(io, μ)
    end
    return nothing
end


###############################################################################
# I <: Tuple

struct TupleProductMeasure{T} <: AbstractProductMeasure
    pars::T
end

export ⊗
⊗(μs::AbstractMeasure...) = productmeasure(μs)

marginals(d::TupleProductMeasure{T}) where {F, T<:Tuple} = d.pars

function Base.show(io::IO, μ::TupleProductMeasure{T}) where {F,T <: Tuple}
    io = IOContext(io, :compact => true)
    print(io, join(string.(marginals(μ)), " ⊗ "))
end

@inline function logdensity(d::TupleProductMeasure, x::Tuple) where {T<:Tuple}
    mapreduce(logdensity, +, d.pars, x)
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::TupleProductMeasure) where {T}
    rand.(d.pars)
end

###############################################################################
# I <: AbstractArray

marginals(d::ProductMeasure{F,A}) where {F,A<:AbstractArray} = mappedarray(d.f, d.pars)

function logdensity(d::ProductMeasure, x)
    mapreduce(logdensity, +, marginals(d), x)
end

function Base.show(io::IO, d::ProductMeasure{F,A}) where {F,A<:AbstractArray}
    io = IOContext(io, :compact => true)
    print(io, "For(")
    print(io, d.f, ", ")
    print(io, d.pars, ")")
end


###############################################################################
# I <: CartesianIndices

function Base.show(io::IO, d::ProductMeasure{F,S,I}) where {F, S, I<:CartesianIndices}
    io = IOContext(io, :compact => true)
    print(io, "For(")
    print(io, d.f, ", ")
    join(io, size(d.pars), ", ")
    print(io, ")")
end


# function Base.rand(rng::AbstractRNG, ::Type{T}, d::ProductMeasure{F,S,I}) where {T,F,I<:CartesianIndices}

# end

###############################################################################
# I <: Base.Generator


export rand!
using Random: rand!, GLOBAL_RNG, AbstractRNG


function logdensity(d::ProductMeasure{F,S,I}, x) where {F, S, I<:Base.Generator}
    sum((logdensity(dj, xj) for (dj, xj) in zip(marginals(d), x)))
end


function Base.rand(rng::AbstractRNG, ::Type{T}, d::ProductMeasure{F,S,I}) where {T,F,S,I<:Base.Generator}
    mar = marginals(d)
    elT = typeof(rand(rng, T, first(mar)))

    sz = size(mar)
    r = ResettableRNG(rng, rand(rng, UInt))
    Base.Generator(s -> rand(r, d.pars.f(s)), d.pars.iter)
end


@propagate_inbounds function Random.rand!(rng::AbstractRNG, d::ProductMeasure, x::AbstractArray)
    # TODO: Generalize this
    T = Float64
    for(j,m) in zip(eachindex(x), marginals(d))
        @inbounds x[j] = rand(rng, T, m)
    end
    return x
end

export rand!
using Random: rand!, GLOBAL_RNG, AbstractRNG

function Base.rand(rng::AbstractRNG, ::Type{T}, d::ProductMeasure) where {T}
    d1 = d.f(first(d.pars))
    rand(rng, T, d, d1)
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::ProductMeasure, d1::AbstractMeasure) where {T}
    mar = marginals(d)
    elT = typeof(rand(rng, T, first(mar)))

    sz = size(mar)
    x = Array{elT, length(sz)}(undef, sz)
    rand!(rng, d, x)
end

# TODO: 
# function Base.rand(rng::AbstractRNG, d::ProductMeasure)
#     return rand(rng, sampletype(d), d)
# end

# function Base.rand(T::Type, d::ProductMeasure)
#     return rand(Random.GLOBAL_RNG, T, d)
# end

# function Base.rand(d::ProductMeasure)
#     T = sampletype(d)
#     return rand(Random.GLOBAL_RNG, T, d)
# end

function sampletype(d::ProductMeasure{A}) where {T,N,A <: AbstractArray{T,N}}
    S = @inbounds sampletype(marginals(d)[1])
    Array{S, N}
end

function sampletype(d::ProductMeasure{<: Tuple}) 
    Tuple{sampletype.(marginals(d))...}
end


# function logdensity(μ::ProductMeasure{Aμ}, x::Ax) where {Aμ <: MappedArray, Ax <: AbstractArray}
#     μ.data
# end

function ConstructionBase.constructorof(::Type{P}) where {F,S,I,P <: ProductMeasure{F,S,I}}
    p -> productmeasure(d.f, p)
end

# function Accessors.set(d::ProductMeasure{N}, ::typeof(params), p) where {N}
#     setproperties(d, NamedTuple{N}(p...))
# end


# function Accessors.set(d::ProductMeasure{F,T}, ::typeof(params), p::Tuple) where {F, T<:Tuple}
#     set.(marginals(d), params, p)
# end

# function logdensity(μ::ProductMeasure, ν::ProductMeasure, x)
#     sum(zip(marginals(μ), marginals(ν), x)) do μ_ν_x
#         logdensity(μ_ν_x...)
#     end
# end

function kernelfactor(μ::ProductMeasure{F,S,<:Fill}) where {F,S}
    k = kernel(first(marginals(μ)))
    (p -> k.f(p)^size(μ), k.ops)
end

function kernelfactor(μ::ProductMeasure{F,S,A}) where {F,S,A<:AbstractArray}
    (p -> set.(marginals(μ), params, p), μ.pars)
end
