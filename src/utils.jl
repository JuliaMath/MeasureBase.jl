const EmptyNamedTuple = NamedTuple{(),Tuple{}}

function Base.show(io::IO, μ::AbstractMeasure)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

showparams(io::IO, ::EmptyNamedTuple) = print(io, "()")
showparams(io::IO, nt::NamedTuple) = print(io, nt)

# @inline function fix(f, x)
#     y = f(x)
#     # Workaround bug https://github.com/JuliaLang/julia/issues/42615
#     while x !== y
#         (x, y) = (y, f(y))
#     end

#     return y
# end

# function constructorof(::Type{T}) where {T} 
#     C = T
#     while C isa UnionAll
#         C = C.body
#     end

#     return C.name.wrapper
# end

constructor(::T) where {T} = constructor(T)

constructor(::Type{T}) where {T} = constructorof(T)

export testvalue
testvalue(μ::AbstractMeasure) = testvalue(basemeasure(μ))

export rootmeasure

basemeasure(μ, x) = basemeasure(μ)

"""
    rootmeasure(μ::AbstractMeasure)

It's sometimes important to be able to find the fix point of a measure under
`basemeasure`. That is, to start with some measure and apply `basemeasure`
repeatedly until there's no change. That's what this does.
"""
@inline function rootmeasure(μ)
    α = basemeasure(μ)
    μ === α && return α
    return rootmeasure(α)
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

@inline function basemeasure_depth(μ::M) where {M<:AbstractMeasure}
    basemeasure_depth(M)
end

@inline basemeasure_depth(::Type{M}) where {M} = static(1) + basemeasure_depth(basemeasure_type(M)) 

export logdensity_tuple

function logdensity_tuple(d, x)
    return (logdensity_def(d, x), basemeasure(d, x), x)
end

function logdensity_tuple(d, (z, x)::MapsTo)
    return (logdensity_def(d, x), basemeasure(d, x), z ↦ x)
end
