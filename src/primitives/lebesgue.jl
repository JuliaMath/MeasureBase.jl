# Lebesgue measure

export Lebesgue

struct LebesgueBase <: PrimitiveMeasure end

@inline mspace_elsize(::LebesgueBase) = ()
@inline mspace_flatsize(::LebesgueBase) = ()
@inline mspace_flatsize(::Type{LebesgueBase}) = ()

massof(::LebesgueBase, s::Interval) = width(s)

testvalue(::LebesgueBase) = 0.0

insupportof(::LebesgueBase, x) = true

insupportof(::LebesgueBase) = Returns(true)

logdensity_rel_def(::LebesgueBase, ::CountingBase, x) = _neg_inf_logd(x)

logdensity_rel_def(::CountingBase, ::LebesgueBase, x) = -_neg_inf_logd(x)

@inline getdof(::LebesgueBase) = static(1)

@inline checked_arg(::LebesgueBase, x::Number) = x

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

@inline mspace_elsize(μ::Lebesgue) = _valueset_elsize(μ.support)
@inline mspace_flatsize(μ::Lebesgue) = _valueset_flatsize(μ.support)
@inline mspace_flatsize(::Type{<:Lebesgue{RealValues}}) = ()
@inline mspace_flatsize(::Type{<:Lebesgue{<:IntervalSets.AbstractInterval}}) = ()

testvalue(::Type{T}, d::Lebesgue) where {T} = testvalue(T, d.support)::T

proxy(d::Lebesgue) = restrict(in(d.support), LebesgueBase())
proxy(::Lebesgue{MeasureBase.RealValues}) = LebesgueBase()

@useproxy Lebesgue

Base.:∘(::typeof(basemeasure), ::Type{Lebesgue}) = LebesgueBase()

Base.show(io::IO, d::Lebesgue) = print(io, "Lebesgue(", d.support, ")")

insupportof(::Lebesgue{RealValues}, ::Real) = true

@inline function logdensityof_impl(μ::Lebesgue, x::Number)
    R = float(typeof(x))
    _checksupport(insupportof(μ, x), zero(R))
end

@inline logdensityof_impl(μ::Lebesgue, x) = _checksupport(insupportof(μ, x), 0.0)

massof(μ::Lebesgue{<:AbstractInterval}) = width(μ.support)

massof(::Lebesgue{RealValues}, s::Interval) = width(s)

# Example: 
# julia> Lebesgue(𝕀)(0.2..5)
# 0.8
function massof(μ::Lebesgue{<:AbstractInterval}, s::Interval)
    a, b = endpoints(μ.support)
    left = max(s.left, a)
    right = min(s.right, b)
    w = right - left
    max(w, zero(w))
end

function smf(μ::Lebesgue{<:AbstractInterval}, x)
    a, b = endpoints(μ.support)
    clamp(x, a, b)
end

invsmf(μ::Lebesgue{<:AbstractInterval}, p) = p

smf(::Lebesgue{<:RealValues}, x) = x
smf(::Lebesgue{<:RealValues}) = identity
invsmf(::Lebesgue{<:RealValues}, x) = x
invsmf(::Lebesgue{<:RealValues}) = identity

smf(::LebesgueBase, x) = x
smf(::LebesgueBase) = identity
invsmf(::LebesgueBase, x) = x
invsmf(::LebesgueBase) = identity


# Coordinate Lebesgue measure on the standard simplex, its own base measure:

const _SimplexLebesgue = Lebesgue{<:StdSimplex}

basemeasure(μ::_SimplexLebesgue) = μ
@inline basemeasure_depth(::_SimplexLebesgue) = static(0)
rootmeasure(μ::_SimplexLebesgue) = μ

@inline getdof(μ::_SimplexLebesgue) = μ.support.n - static(1)
@inline fast_dof(μ::_SimplexLebesgue) = getdof(μ)

massof(μ::_SimplexLebesgue) = inv(prod(float, 1:(dynamic(μ.support.n) - 1); init = 1.0))

gentype(::_SimplexLebesgue) = Vector{Float64}

testvalue(::Type{T}, μ::_SimplexLebesgue) where {T} = testvalue(T, μ.support)

insupportof(μ::_SimplexLebesgue) = Base.Fix1(insupportof, μ)
localmeasure(μ::_SimplexLebesgue, ::Any) = μ
transportmeasure(μ::_SimplexLebesgue, ::Any) = μ

@inline mspace_ndims(::Type{<:_SimplexLebesgue}) = 1

@propagate_inbounds function checked_arg(μ::_SimplexLebesgue, x::AbstractVector)
    @boundscheck length(x) == μ.support.n || _throw_size_mismatch()
    return x
end
checked_arg(::_SimplexLebesgue, ::Any) = _throw_size_mismatch()

logdensity_def(::_SimplexLebesgue, x) = zero(_logd_numtype(x))
unsafe_logdensityof(::_SimplexLebesgue, x) = zero(_logd_numtype(x))
logdensityof_impl(μ::_SimplexLebesgue, x::AbstractVector) = batched_logdensityof_impl(μ, x)

function batched_logdensityof_impl(::_SimplexLebesgue, X::AbstractArray)
    T = _logd_numtype(X)
    ifelse.(_in_stdsimplex(X), zero(T), T(-Inf))
end
batched_logdensity_def(::_SimplexLebesgue, X::AbstractArray) = _zero_logd_batch(X, static(1))

# The normalized law is the uniform distribution on the simplex, the
# normalized vector of independent standard exponential variates:
function batched_rand_impl(ctx::GenContext, μ::_SimplexLebesgue, sz::SizeLike)
    E = _randexp_bulk(ctx, (μ.support.n, size_tuple(sz)...))
    E ./ sum(E; dims = 1)
end
