"""A `NamedTuple` with no element."""
const EmptyNamedTuple = NamedTuple{(),Tuple{}}

function Base.show(io::IO, μ::AbstractMeasure)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

showparams(io::IO, ::EmptyNamedTuple) = print(io, "()")
showparams(io::IO, nt::NamedTuple) = print(io, nt)

"""
    fix(f, x)

Find a fixed point of `f` starting from point `x` by iterative application.
"""
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

"""
    testvalue(μ::AbstractMeasure)

Return a typical value from the domain of `μ`.
"""
testvalue(μ::AbstractMeasure) = testvalue(basemeasure(μ))

"""
    rootmeasure(μ::AbstractMeasure)

Find the fixed point of `basemeasure` starting from measure `μ`.
"""
rootmeasure(μ::AbstractMeasure) = fix(basemeasure, μ)

# Base on the Tricks.jl README

struct Iterable end

struct NonIterable end

"""
    isiterable(::Type{T})

Determine whether type `T` is an iterable by checking for an `iterate` method.
"""
function isiterable(::Type{T}) where {T}
    return static_hasmethod(iterate, Tuple{T}) ? Iterable() : NonIterable()
end

functioninstance(::Type{F}) where {F<:Function} = F.instance

# See https://github.com/cscherrer/KeywordCalls.jl/issues/22
@inline instance_type(f::F) where {F<:Function} = F
@inline instance_type(f::UnionAll) = Type{f}
