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
    inbounds = isnonnegative
    constℓ = static(logtwo)
    varℓ = Returns(0.0)
    base = basemeasure(unhalf(μ))
    FactoredBase(inbounds, constℓ, varℓ, base)
end

function basemeasure_type(::Type{Half{M}}) where {M}
    B = basemeasure_type(M)
    FactoredBase{typeof(isnonnegative), StaticFloat64{logtwo}, Returns{Float64}, B}
end

function Base.rand(rng::AbstractRNG, T::Type, μ::Half)
    return abs(rand(rng, T, unhalf(μ)))
end

logdensity_def(μ::Half, x) = logdensity_def(unhalf(μ), x)
