export ParameterizedMeasure

abstract type ParameterizedMeasure{N} <: AbstractMeasure end

function Base.getproperty(μ::ParameterizedMeasure{N}, prop::Symbol) where {N}
    return getproperty(getfield(μ, :par), prop)
end

function Base.propertynames(μ::ParameterizedMeasure{N}) where {N}
    return N
end

function Base.show(io::IO, μ::ParameterizedMeasure{()})
    io = IOContext(io, :compact => true)
    print(io, nameof(typeof(μ)), "()")
end

function Base.show(io::IO, μ::ParameterizedMeasure{N}) where {N}
    io = IOContext(io, :compact => true)
    print(io, nameof(typeof(μ)))
    print(io, getfield(μ, :par))
end

# Allow things like
#
# julia> Normal{(:μ,)}(2)
# Normal(μ = 2,)
#
export kleisli

function kleisli(::Type{P}) where {N,P<:ParameterizedMeasure{N}}
    function(args...) P(args...) end
end

function (::Type{P})(args...) where {N,P<:ParameterizedMeasure{N}}
    C = constructorof(P)
    return C(NamedTuple{N}(args...))
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
# kleislifactor

function kleislifactor(::Type{P}) where {N,P<:ParameterizedMeasure{N}}
    (constructorof(P), N)
end

function kleislifactor(::P) where {N,P<:ParameterizedMeasure{N}}
    (constructorof(P), N)
end
