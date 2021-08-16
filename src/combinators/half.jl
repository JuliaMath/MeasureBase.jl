export Half

struct Half{M} <: AbstractMeasure
    parent::M
end

function Base.show(io::IO, μ::Half)
    print(io, "Half")
    show(io, μ.parent)
end

unhalf(μ::Half) = μ.parent

basemeasure(μ::Half) = basemeasure(μ.parent)

function Base.rand(rng::AbstractRNG, T::Type, μ::Half)
    return abs(rand(rng, T, unhalf(μ)))
end
