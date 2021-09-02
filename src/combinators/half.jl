export Half

struct Half{M} <: AbstractMeasure
    parent::M
end

function Base.show(io::IO, μ::Half)
    print(io, "Half")
    show(io, μ.parent)
end

unhalf(μ::Half) = μ.parent

basemeasure(μ::Half) = WeightedMeasure(logtwo, basemeasure(unhalf(μ)))

function Base.rand(rng::AbstractRNG, T::Type, μ::Half)
    return abs(rand(rng, T, unhalf(μ)))
end

logdensity(μ::Half, x) = x > 0 ? logdensity(unhalf(μ), x) : -Inf
