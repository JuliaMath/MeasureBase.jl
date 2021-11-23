export Half

struct Half{M} <: AbstractMeasure
    parent::M
end

function Base.show(io::IO, μ::Half)
    print(io, "Half")
    show(io, μ.parent)
end

unhalf(μ::Half) = μ.parent

@inline function basemeasure(μ::Half)
    inbounds(x) = x > 0
    constℓ = logtwo
    varℓ() = 0.0
    base = basemeasure(unhalf(μ))
    FactoredBase(inbounds, constℓ, varℓ, base)
end

function Base.rand(rng::AbstractRNG, T::Type, μ::Half)
    return abs(rand(rng, T, unhalf(μ)))
end

logdensity_def(μ::Half, x) = logdensity_def(unhalf(μ), x)
