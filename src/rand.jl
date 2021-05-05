Base.rand(d::AbstractMeasure) = rand(Random.GLOBAL_RNG, Float64, d)

Base.rand(T::Type, μ::AbstractMeasure) = rand(Random.GLOBAL_RNG, T, μ)

Base.rand(rng::AbstractRNG, d::AbstractMeasure) = rand(rng, Float64, d)

# Base.rand(rng::AbstractRNG, T::Type, d::ParameterizedMeasure) = rand(rng, distproxy(d))

# TODO: Make this easily configurable, e.g. we should be able to use `collect` or `tuple` at least
# Maybe RandomExtensions can help here?
function Base.rand(rng::AbstractRNG, T::Type, d::ProductMeasure)
    collect((rand(rng, T, dn) for dn in d.data))
end


@inline Random.rand!(d::AbstractMeasure, arr::AbstractArray) = rand!(GLOBAL_RNG, d, arr)
