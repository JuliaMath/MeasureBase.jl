
export Dirac

struct Dirac{X} <: PrimitiveMeasure
    x::X
end

gentype(μ::Dirac{X}) where {X} = X

function (μ::Dirac{X})(s) where {X}
    μ.x ∈ s && return 1
    return 0
end

density_def(μ::Dirac, x) = x == μ.x

logdensity_def(μ::Dirac, x) = (x == μ.x) ? 0.0 : -Inf

Base.rand(::Random.AbstractRNG, T::Type, μ::Dirac) = μ.x

export dirac

dirac(d::AbstractMeasure) = Dirac(rand(d))

@inline function logdensity_def(μ::Dirac{M}, ν::Dirac{M}, x) where {M}
    logdensity_def(μ, x) - logdensity_def(ν, x)
end

testvalue(d::Dirac) = d.x
