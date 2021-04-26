
export ParameterizedMeasure
abstract type ParameterizedMeasure{N} <: AbstractMeasure end

function Base.getproperty(μ::ParameterizedMeasure{N}, prop::Symbol) where {N}
    return getproperty(getfield(μ, :par), prop)
end

function Base.propertynames(μ::ParameterizedMeasure{N}) where {N}
    return N
end

function Base.show(io::IO, μ::ParameterizedMeasure{()}) 
    print(io, nameof(typeof(μ)), "()")
end

function Base.show(io::IO, μ::ParameterizedMeasure{N}) where {N}
    io = IOContext(io, :compact => true)
    print(io, nameof(typeof(μ)))
    print(io, getfield(μ,:par))
end

# e.g. Normal(;μ=μ, σ=σ) = Normal((μ=μ, σ=σ))
(M::Type{<: ParameterizedMeasure})(; kwargs...) = M(paramsort((; kwargs...)))

(M::Type{<: ParameterizedMeasure})(::Tuple{}) = M(NamedTuple())

function (M::Type{P})(args...) where {N, P <: ParameterizedMeasure{N}}
    constructor = M.body.name.wrapper
    return constructor(NamedTuple{N}(args...))
end

export params

params(μ::ParameterizedMeasure) = getfield(μ, :par)
