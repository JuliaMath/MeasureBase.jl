
export Dirac

struct Dirac{X} <: AbstractMeasure
    x::X
end

function Pretty.tile(d::Dirac)
    Pretty.literal("Dirac(") * Pretty.tile(d.x) * Pretty.literal(")")
end

gentype(μ::Dirac{X}) where {X} = X

function (μ::Dirac{X})(s) where {X}
    μ.x ∈ s && return 1
    return 0
end

basemeasure(d::Dirac) = CountingMeasure()

logdensity_def(μ::Dirac, x) = 0.0

Base.rand(::Random.AbstractRNG, T::Type, μ::Dirac) = μ.x

export dirac

dirac(d::AbstractMeasure) = Dirac(rand(d))

testvalue(d::Dirac) = d.x

insupport(d::Dirac, x) = x == d.x

@inline getdof(::Dirac) = static(0)

@propagate_inbounds function checked_var(μ::Dirac, x)
    @boundscheck insupport(μ, x) || throw(ArgumentError("Invalid variate for measure"))
    x
end

@inline vartransform_def(ν::Dirac, ::PowerMeasure{<:MeasureBase.StdMeasure}, ::Any) = ν.x

@inline function vartransform_def(ν::PowerMeasure{<:MeasureBase.StdMeasure}, ::Dirac, ::Any)
    Zeros{Bool}(map(_ -> 0, ν.axes))
end
