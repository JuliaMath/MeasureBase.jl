export ParameterizedMeasure

"""
    abstract type ParameterizedMeasure{N} <: AbstractMeasure

Supertype of measures with named parameters `N`, stored as a `NamedTuple`
in a field `par`.

Subtypes `M{N,T} <: ParameterizedMeasure{N}` can be constructed as `M(nt)`
or with positional parameters `M{N}(args...)`, and support `paramsof`,
[`paramnames`](@ref) and property access to the parameters.
"""
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
    result *= Pretty.literal(sprint(show, par; context = :compact => true))
    result
end

function Pretty.tile(d::ParameterizedMeasure{()})
    result = Pretty.literal(nameof(typeof(d)))
    par = getfield(d, :par)
    result *= Pretty.literal("()")
    result
end

function (::Type{P})(nt::NamedTuple{K,T}) where {K,T,N,P<:ParameterizedMeasure{N}}
    C = constructorof(P)
    arg = NamedTuple{N}(nt)
    return _parameterized(C, arg)
end

function _parameterized(::Type{C}, arg::NamedTuple{N,T}) where {C,N,T}
    return C(arg)::C{N,T}
end

function (::Type{P})(args...) where {N,P<:ParameterizedMeasure{N}}
    C = constructorof(P)
    return C(NamedTuple{N}(args))::C{N,typeof(args)}
end

function ConstructionBase.setproperties(
    d::P,
    nt::NamedTuple,
) where {P<:ParameterizedMeasure}
    return constructorof(P)(merge(paramsof(d), nt))
end

###############################################################################
# paramsof

export paramsof

"""
`paramsof(μ)` returns the parameters of a measure `μ`, as a `NamedTuple`. The
default method is
```
paramsof(μ) = NamedTuple()
```

See also `paramnames`
"""
function paramsof end

@deprecate params(μ) paramsof(μ)
@deprecate params(μ, constraints) paramsof(μ, constraints)

paramsof(μ::ParameterizedMeasure) = getfield(μ, :par)

function paramsof(μ::AbstractMeasure, constraints::NamedTuple{C}) where {C}
    NamedTuple{paramnames(μ, constraints)}(paramsof(μ))
end

paramsof(μ) = NamedTuple()

###############################################################################
# paramnames

export paramnames

"""
`paramnames(μ)` returns the names of the parameters of `μ`. This is equivalent to 
```
paramnames(μ) == (keys ∘ paramsof)(μ)
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

See also `paramsof`
"""
function paramnames end

paramnames(μ::M) where {M} = paramnames(M)

paramnames(::Type{PM}) where {N,PM<:ParameterizedMeasure{N}} = N

paramnames(μ::AbstractMeasure) = propertynames(μ)

paramsof(::Type{PM}) where {N,PM<:ParameterizedMeasure{N}} = N

function paramnames(μ, constraints::NamedTuple{N}) where {N}
    tuple((k for k in paramnames(μ) if k ∉ N)...)
end
