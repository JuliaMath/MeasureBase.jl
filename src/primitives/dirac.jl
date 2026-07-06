
export Dirac

struct Dirac{X} <: AbstractMeasure
    x::X
end

function Pretty.tile(d::Dirac)
    Pretty.literal("Dirac(") * Pretty.tile(d.x) * Pretty.literal(")")
end

Base.:(==)(a::Dirac, b::Dirac) = a.x == b.x
Base.isapprox(a::Dirac, b::Dirac; kwargs...) = isapprox(a.x, b.x; kwargs...)

gentype(μ::Dirac{X}) where {X} = X

function (μ::Dirac{X})(s) where {X}
    μ.x ∈ s && return 1
    return 0
end

basemeasure(d::Dirac) = CountingBase()

massof(::Dirac) = static(1.0)

function logdensityof_impl(μ::Dirac, x::Real)
    R = float(typeof(x))
    insupport(μ, x) ? zero(R) : R(-Inf)
end

logdensityof_impl(μ::Dirac, x) = insupport(μ, x) ? 0.0 : -Inf

logdensity_def(::Dirac, x::Real) = zero(float(typeof(x)))
logdensity_def(::Dirac, x) = 0.0

Base.rand(::Random.AbstractRNG, T::Type, μ::Dirac) = μ.x

export dirac

dirac(d::AbstractMeasure) = Dirac(rand(d))

insupport(d::Dirac, x) = x == d.x

@inline getdof(::Dirac) = static(0)

@inline mspace_elsize(μ::Dirac) = maybestatic_size(μ.x)

@propagate_inbounds function checked_arg(μ::Dirac, x)
    @boundscheck insupport(μ, x) || throw(ArgumentError("Invalid variate for measure"))
    x
end
