const EmptyNamedTuple = NamedTuple{(),Tuple{}}

"""
    function proxy end

It's often useful to delegate methods like `logdensity` and `basemeasure` to
those of a different measure. For example, a `Normal{(:μ,:σ)}` is equivalent to
an affine transformation of a `Normal{()}`.

We _could_ just have calls like `Normal(μ=2,σ=4)` directly construct a
transformed measure, but this would make dispatch awkward.
"""
proxy(μ) = μ

function Base.show(io::IO, μ::AbstractMeasure)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

showparams(io::IO, ::EmptyNamedTuple) = print(io, "()")
showparams(io::IO, nt::NamedTuple) = print(io, nt)

@inline function fix(f, x)
    y = f(x)
    # Workaround bug https://github.com/JuliaLang/julia/issues/42615
    while x !== y
        (x, y) = (y, f(y))
    end

    return y
end

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

using MLStyle

export basemeasure_depth

basemeasure_depth(μ::M) where {M<:AbstractMeasure} = basemeasure_depth(M)

@inline function basemeasure_depth(μ::M) where {M}
    static(1) + basemeasure_depth(basemeasure(μ))
end

@inline basemeasure_depth(::PrimitiveMeasure) = static(0)

export logdensity_tuple

logdensity_tuple(d, x) = logdensity_tuple(MeasureBase.proxy(d), x)

function logdensity_tuple(d, (z,x)::MapsTo)
    prox = MeasureBase.proxy(d)
    return (logdensity(prox, x), basemeasure(prox, x), z ↦ x)
end