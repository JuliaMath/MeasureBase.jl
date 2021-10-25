
export Dirac

struct Dirac{X} <: PrimitiveMeasure
    x::X
end

function Pretty.quoteof(d::Dirac)
    q = Pretty.quoteof(d.x)
    :(Dirac($q))
end

sampletype(μ::Dirac{X}) where {X} = X

function (μ::Dirac{X})(s) where {X}
    μ.x ∈ s && return 1
    return 0
end

density(μ::Dirac, x) = x == μ.x

logdensity(μ::Dirac, x) = (x == μ.x) ? 0.0 : -Inf


Base.rand(::Random.AbstractRNG, T::Type, μ::Dirac) = μ.x


export dirac

dirac(d::AbstractMeasure) = Dirac(rand(d))

function logdensity(μ::Dirac{M}, ν::Dirac{M}, x) where {M} 
    logdensity(μ,x) - logdensity(ν,x)
end

testvalue(d::Dirac) = d.x
