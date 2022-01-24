export Half

struct Half{M} <: AbstractMeasure
    parent::M
end

function Base.show(io::IO, μ::Half)
    print(io, "Half")
    show(io, μ.parent)
end

unhalf(μ::Half) = μ.parent

isnonnegative(x) = x ≥ 0.0

@inline function basemeasure(μ::Half)
    constℓ = static(logtwo)
    varℓ = Returns(static(0.0))
    base = basemeasure(unhalf(μ))
    FactoredBase(constℓ, varℓ, base)
end

function Base.rand(rng::AbstractRNG, T::Type, μ::Half)
    return abs(rand(rng, T, unhalf(μ)))
end

logdensity_def(μ::Half, x) = logdensity_def(unhalf(μ), x)

@inline function insupport(d::Half, x)
    ifelse(isnonnegative(x), insupport(unhalf(d), x), false)
end