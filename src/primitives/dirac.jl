
export Dirac

struct Dirac{X} <: AbstractMeasure
    x::X
end

gentype(μ::Dirac{X}) where {X} = X

function (μ::Dirac{X})(s) where {X}
    μ.x ∈ s && return 1
    return 0
end

basemeasure(d::Dirac) = CountingMeasure()

basemeasure_depth(::Dirac) = static(1)
basemeasure_depth(::Type{D}) where {D<:Dirac} = static(1)



density_def(μ::Dirac, x) = x == μ.x

logdensity_def(μ::Dirac, x) = (x == μ.x) ? 0.0 : -Inf

Base.rand(::Random.AbstractRNG, T::Type, μ::Dirac) = μ.x

export dirac

dirac(d::AbstractMeasure) = Dirac(rand(d))

testvalue(d::Dirac) = d.x
