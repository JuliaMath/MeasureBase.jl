export ParameterizedMeasure

abstract type ParameterizedMeasure{N} <: AbstractMeasure end

function Base.getproperty(μ::ParameterizedMeasure{N}, prop::Symbol) where {N}
    return getproperty(getfield(μ, :par), prop)
end

function Base.propertynames(μ::ParameterizedMeasure{N}) where {N}
    return N
end

function Pretty.tile(d::ParameterizedMeasure)
    result = Pretty.literal(nameof(typeof(d)))
    par = getfield(d, :par)
    result *= Pretty.literal(sprint(show, par; context=:compact => true))
    result
end

function Pretty.tile(d::ParameterizedMeasure{()})
    result = Pretty.literal(nameof(typeof(d)))
    par = getfield(d, :par)
    result *= Pretty.literal("()")
    result
end

# function Base.show(io::IO, μ::ParameterizedMeasure{()})
#     io = IOContext(io, :compact => true)
#     print(io, nameof(typeof(μ)), "()")
# end

# function Base.show(io::IO, μ::ParameterizedMeasure{N}) where {N}
#     io = IOContext(io, :compact => true)
#     print(io, nameof(typeof(μ)))
#     print(io, getfield(μ, :par))
# end

# Allow things like
#
# julia> Normal{(:μ,)}(2)
# Normal(μ = 2,)
#

function kernel(::Type{P}) where {N,P<:ParameterizedMeasure{N}}
    C = constructorof(P)
    _kernel(C, Val(N))
end

@inline function _kernel(::Type{C}, ::Val{N}) where {C,N}
    @inline function f(args::T) where {T<:Tuple}
        C(NamedTuple{N,T}(args))::C{N,T}
    end

    @inline function f(arg::T) where {T}
        C(NamedTuple{N,Tuple{T}}((arg,)))::C{N,Tuple{T}}
    end

    f
end

function (::Type{P})(nt::NamedTuple) where {N,P<:ParameterizedMeasure{N}}
    C = constructorof(P)
    arg = NamedTuple{N}(nt)
    return C(arg)::P
end

function (::Type{P})(args...) where {N,P<:ParameterizedMeasure{N}}
    C = constructorof(P)
    return C(NamedTuple{N}(args))::C{N,typeof(args)}
end

(::Type{P})(; kwargs...) where {P<:ParameterizedMeasure} = P(NamedTuple(kwargs))

function ConstructionBase.setproperties(
    d::P,
    nt::NamedTuple,
) where {P<:ParameterizedMeasure}
    return constructorof(P)(merge(params(d), nt))
end

###############################################################################
# params

export params

"""
`params(μ)` returns the parameters of a measure `μ`, as a `NamedTuple`. The
default method is
```
params(μ) = NamedTuple()
```

See also `paramnames`
"""
function params end

params(μ::ParameterizedMeasure) = getfield(μ, :par)

function params(μ::AbstractMeasure, constraints::NamedTuple{C}) where {C}
    NamedTuple{paramnames(μ, constraints)}(params(μ))
end

params(μ) = NamedTuple()

###############################################################################
# paramnames

export paramnames

"""
`paramnames(μ)` returns the names of the parameters of `μ`. This is equivalent to 
```
paramnames(μ) == (keys ∘ params)(μ)
```
but depends only on the type. In particular, the default implementation is
```
paramnames(μ::M) where {M} = paramnames(M)
```

New `ParameterizedMeasure`s will automatically have a `paramnames` method. For
other measures, this method is optional, but can be added by defining
```
paramnames(::Type{M}) where {M} = ...
```

See also `params`
"""
function paramnames end

paramnames(μ::M) where {M} = paramnames(M)

paramnames(::Type{PM}) where {N,PM<:ParameterizedMeasure{N}} = N

paramnames(μ::AbstractMeasure) = propertynames(μ)

params(::Type{PM}) where {N,PM<:ParameterizedMeasure{N}} = N

function paramnames(μ, constraints::NamedTuple{N}) where {N}
    tuple((k for k in paramnames(μ) if k ∉ N)...)
end

###############################################################################
# kernelfactor

function kernelfactor(::Type{P}) where {N,P<:ParameterizedMeasure{N}}
    (constructorof(P), N)
end

function kernelfactor(::P) where {N,P<:ParameterizedMeasure{N}}
    (constructorof(P), N)
end
