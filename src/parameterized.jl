"""Abstract supertype for measures defined by a parameter `par::NamedTuple{N}`."""
abstract type ParameterizedMeasure{N} <: AbstractMeasure end

function Base.getproperty(μ::ParameterizedMeasure{N}, prop::Symbol) where {N}
    return getproperty(getfield(μ, :par), prop)
end

function Base.propertynames(::ParameterizedMeasure{N}) where {N}
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
function (::Type{P})(args...) where {N,P<:ParameterizedMeasure{N}}
    C = constructorof(P)
    return C(NamedTuple{N}(args...))
end

(::Type{P})(; kwargs...) where {P<:ParameterizedMeasure} = P(NamedTuple(kwargs))

function ConstructionBase.setproperties(
    d::P,
    nt::NamedTuple
) where {P<:ParameterizedMeasure}
    return constructorof(P)(merge(params(d), nt))
end

###############################################################################
# params

"""
    params(μ::AbstractMeasure)

Get the parameters of measure `μ` as a `NamedTuple`.
"""
params(::AbstractMeasure) = NamedTuple()
params(μ::ParameterizedMeasure) = getfield(μ, :par)
params(::Type{PM}) where {N,PM<:ParameterizedMeasure{N}} = N

function params(μ::AbstractMeasure, constraints::NamedTuple{C}) where {C}
    NamedTuple{paramnames(μ, constraints)}(params(μ))
end

###############################################################################
# paramnames

paramnames(μ) = paramnames(typeof(μ))
paramnames(::Type{PM}) where {N,PM<:ParameterizedMeasure{N}} = N

"""
    paramnames(μ::AbstractMeasure)

Get the names of the parameters for measure `μ`.
"""
paramnames(μ::AbstractMeasure) = propertynames(μ)

function paramnames(μ, constraints::NamedTuple{N}) where {N}
    tuple((k for k in paramnames(μ) if !(k in N))...)
end

###############################################################################
# kernelfactor

function kernelfactor(::Type{P}) where {N,P<:ParameterizedMeasure{N}}
    (constructorof(P), N)
end

function kernelfactor(::P) where {N,P<:ParameterizedMeasure{N}}
    (constructorof(P), N)
end
