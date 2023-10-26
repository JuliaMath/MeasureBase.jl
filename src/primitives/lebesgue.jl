# Lebesgue measure

export Lebesgue

struct LebesgueBase <: PrimitiveMeasure end

massof(::LebesgueBase, s::Interval) = width(s)

testvalue(::LebesgueBase) = 0.0

insupport(::LebesgueBase, x) = true

insupport(::LebesgueBase) = Returns(true)

logdensity_def(::LebesgueBase, ::CountingBase, x) = -Inf

logdensity_def(::CountingBase, ::LebesgueBase, x) = Inf

@inline getdof(::LebesgueBase) = static(1)

@inline checked_arg(::LebesgueBase, x::Real) = x

@propagate_inbounds function checked_arg(::LebesgueBase, x::Any)
    @boundscheck throw(ArgumentError("Invalid variate type for measure"))
end

massof(::LebesgueBase) = static(Inf)

function _default_massof_impl(m, s::AbstractInterval, ::LebesgueBase)
    mass = massof(m)
    nu = mass * StdUniform()
    f = transport_to(nu, m)
    a = f(leftendpoint(s))
    b = f(rightendpoint(s))
    return mass * abs(b - a)
end

##########################################################
struct Lebesgue{T} <: AbstractMeasure
    support::T
end

function Pretty.tile(μ::Lebesgue)
    Pretty.list_layout([Pretty.tile(μ.support)]; prefix = :Lebesgue)
end

gentype(::Lebesgue) = Float64

Lebesgue() = Lebesgue(ℝ)

testvalue(::Type{T}, d::Lebesgue) where {T} = testvalue(T, d.support)::T

proxy(d::Lebesgue) = restrict(in(d.support), LebesgueBase())
proxy(::Lebesgue{MeasureBase.RealNumbers}) = LebesgueBase()

@useproxy Lebesgue

Base.:∘(::typeof(basemeasure), ::Type{Lebesgue}) = LebesgueBase()

Base.show(io::IO, d::Lebesgue) = print(io, "Lebesgue(", d.support, ")")

insupport(μ::Lebesgue, x) = x ∈ μ.support

insupport(::Lebesgue{RealNumbers}, ::Real) = true

massof(::Lebesgue{RealNumbers}, s::Interval) = width(s)

# Example: 
# julia> Lebesgue(𝕀)(0.2..5)
# 0.8
function massof(μ::Lebesgue{<:BoundedReals}, s::Interval)
    a = μ.support.lower
    b = μ.support.upper
    left = max(s.left, a)
    right = min(s.right, b)
    w = right - left
    max(w, zero(w))
end

function smf(μ::Lebesgue{<:BoundedReals}, x)
    clamp(x, μ.support.lower, μ.support.upper)
end

smf(::Lebesgue{RealNumbers}, x) = x
smf(::Lebesgue{RealNumbers}) = identity
invsmf(::Lebesgue{RealNumbers}, x) = x
invsmf(::Lebesgue{RealNumbers}) = identity

smf(::LebesgueBase, x) = x
smf(::LebesgueBase) = identity
invsmf(::LebesgueBase, x) = x
invsmf(::LebesgueBase) = identity
