import ConstructionBase

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
    print(io, getfield(μ,:par))
end


# Allow things like
#
# julia> Normal{(:μ,)}(2)
# Normal(μ = 2,)
#
function (::Type{P})(args...) where {N, P <: ParameterizedMeasure{N}}
    C = constructorof(P)
    return C(NamedTuple{N}(args...))
end

(::Type{P})(;kwargs...) where {P <: ParameterizedMeasure} = P(NamedTuple(kwargs))

export params

params(::Type{PM}) where {N, PM<:ParameterizedMeasure{N}} = N

function params(::Type{M}, constraints::NamedTuple{N2}) where {N1, N2, M<: ParameterizedMeasure{N1}} 
    tuple((k for k in N1 if k ∉ N2)...)
end


params(μ) = ()


function ConstructionBase.setproperties(d::P, nt::NamedTuple) where {P<:ParameterizedMeasure}
    return constructorof(P)(merge(params(d), nt)) 
end
