
export Dirac

struct Dirac{X} <: AbstractMeasure
    x::X
end

function Pretty.quoteof(d::Dirac)
    q = Pretty.quoteof(d.x)
    :(Dirac($q))
end

gentype(μ::Dirac{X}) where {X} = X

function (μ::Dirac{X})(s) where {X}
    μ.x ∈ s && return 1
    return 0
end

basemeasure(d::Dirac) = CountingMeasure()

tbasemeasure_type(::Type{Dirac}) = CountingMeasure

density_def(μ::Dirac, x) = x == μ.x

logdensity_def(μ::Dirac, x) = (x == μ.x) ? 0.0 : -Inf

Base.rand(::Random.AbstractRNG, T::Type, μ::Dirac) = μ.x

export dirac

dirac(d::AbstractMeasure) = Dirac(rand(d))

testvalue(d::Dirac) = d.x

insupport(d::Dirac, x) = x == d.x
