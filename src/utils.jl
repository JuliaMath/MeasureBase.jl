const EmptyNamedTuple = NamedTuple{(),Tuple{}}

function Base.show(io::IO, μ::AbstractMeasure)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

showparams(io::IO, ::EmptyNamedTuple) = print(io, "()")
showparams(io::IO, nt::NamedTuple) = print(io, nt)


export testvalue
testvalue(μ::AbstractMeasure) = testvalue(basemeasure(μ))
testvalue(::Type{T}) where {T} = zero(T)

export rootmeasure

basemeasure(μ, x) = basemeasure(μ)

"""
    rootmeasure(μ::AbstractMeasure)

It's sometimes important to be able to find the fix point of a measure under
`basemeasure`. That is, to start with some measure and apply `basemeasure`
repeatedly until there's no change. That's what this does.
"""
@inline function rootmeasure(μ)
    n = basemeasure_depth(μ)
    _rootmeasure(μ, static(n))
end

@generated function _rootmeasure(μ, ::StaticInt{n}) where {n}
    q = quote end
    foreach(1:n) do _
        push!(q.args, :(μ = basemeasure(μ)))
    end
    return q
end

# Base on the Tricks.jl README
using Tricks
struct Iterable end
struct NonIterable end
isiterable(::Type{T}) where {T} =
    static_hasmethod(iterate, Tuple{T}) ? Iterable() : NonIterable()

functioninstance(::Type{F}) where {F<:Function} = F.instance

# See https://github.com/cscherrer/KeywordCalls.jl/issues/22
@inline instance_type(f::F) where {F<:Function} = F
@inline instance_type(f::UnionAll) = Type{f}

export basemeasure_depth

@inline @constprop :aggressive function basemeasure_type(μ::M) where {M<:AbstractMeasure}
    return tbasemeasure_type(M)
end

@inline @constprop :aggressive function basemeasure_depth(μ::M) where {M<:AbstractMeasure}
    return tbasemeasure_depth(M)
end

@inline @constprop :aggressive tbasemeasure_depth(::Type{M}) where M = tbasemeasure_depth(M, tbasemeasure_type(M), static(0))

@generated function tbasemeasure_depth(::Type{M}, ::Type{B}, S::StaticInt{N}) where {M,B,N}
    M === B && return static(N)
    return tbasemeasure_depth(B, tbasemeasure_type(B), static(N+1))
end