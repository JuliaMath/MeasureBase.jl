# TODO: Compare with ChangesOfVariables.jl

abstract type AbstractTransformedMeasure <: AbstractMeasure end

abstract type AbstractPushforward <: AbstractTransformedMeasure end

abstract type AbstractPullback <: AbstractTransformedMeasure end

function gettransform(::AbstractTransformedMeasure) end

function params(::AbstractTransformedMeasure) end

function paramnames(::AbstractTransformedMeasure) end

function parent(::AbstractTransformedMeasure) end

export PushforwardMeasure

"""
    struct PushforwardMeasure{FF,IF,MU,VC<:TransformVolCorr} <: AbstractPushforward
        f :: FF
        inv_f :: IF
        origin :: MU
        volcorr :: VC
    end
"""
struct PushforwardMeasure{FF,IF,M,VC<:TransformVolCorr} <: AbstractPushforward
    f::FF
    inv_f::IF
    origin::M
    volcorr::VC
end


gettransform(ν::PushforwardMeasure) = ν.f
parent(ν::PushforwardMeasure) = ν.origin

function transport_def(ν::PushforwardMeasure{FF,IF,M}, μ::M, x) where {FF,IF,M}
    if μ == parent(ν)
        return ν.f(x)
    else
        invoke(transport_def, Tuple{Any, Any, Any}, ν, μ, x)
    end
end

function transport_def(μ::M, ν::PushforwardMeasure{FF,IF,M}, y) where {FF,IF,M}
    if μ == parent(ν)
        return ν.inv_f(y)
    else
        invoke(transport_def, Tuple{Any, Any, Any}, μ, ν, y)
    end
end

function Pretty.tile(ν::PushforwardMeasure)
    Pretty.list_layout(Pretty.tile.([ν.f, ν.inv_f, ν.origin]); prefix = :PushforwardMeasure)
end

@inline function logdensity_def(
    ν::PushforwardMeasure{FF,IF,M,<:WithVolCorr},
    y,
) where {FF,IF,M}
    x_orig, inv_ladj = with_logabsdet_jacobian(ν.inv_f, y)
    logd_orig = logdensity_def(ν.origin, x_orig)
    logd = float(logd_orig + inv_ladj)
    neginf = oftype(logd, -Inf)
    return ifelse(
        # Zero density wins against infinite volume:
        (isnan(logd) && logd_orig == -Inf && inv_ladj == +Inf) ||
        # Maybe  also for (logd_orig == -Inf) && isfinite(inv_ladj) ?
        # Return constant -Inf to prevent problems with ForwardDiff:
        (isfinite(logd_orig) && (inv_ladj == -Inf)),
        neginf,
        logd,
    )
end

@inline function logdensity_def(
    ν::PushforwardMeasure{FF,IF,M,<:NoVolCorr},
    y,
) where {FF,IF,M}
    x_orig = to_origin(ν, y)
    return logdensity_def(ν.origin, x_orig)
end

insupport(ν::PushforwardMeasure, y) = insupport(transport_origin(ν), to_origin(ν, y))

testvalue(::Type{T}, ν::PushforwardMeasure) where {T} = from_origin(ν, testvalue(T, transport_origin(ν)))

@inline function basemeasure(ν::PushforwardMeasure)
    PushforwardMeasure(ν.f, ν.inv_f, basemeasure(transport_origin(ν)), NoVolCorr())
end

_pushfwd_dof(::Type{MU}, ::Type, dof) where {MU} = NoDOF{MU}()
_pushfwd_dof(::Type{MU}, ::Type{<:Tuple{Any,Real}}, dof) where {MU} = dof

# Assume that DOF are preserved if with_logabsdet_jacobian is functional:
@inline function getdof(ν::MU) where {MU<:PushforwardMeasure}
    T = Core.Compiler.return_type(testvalue, Tuple{typeof(ν.origin)})
    R = Core.Compiler.return_type(with_logabsdet_jacobian, Tuple{typeof(ν.f),T})
    _pushfwd_dof(MU, R, getdof(ν.origin))
end

# Bypass `checked_arg`, would require potentially costly transformation:
@inline checked_arg(::PushforwardMeasure, x) = x

# @inline transport_origin(ν::PushforwardMeasure) = ν.origin
# @inline from_origin(ν::PushforwardMeasure, x) = ν.f(x)
# @inline to_origin(ν::PushforwardMeasure, y) = ν.inv_f(y)

@inline transport_origin(μ::PushforwardMeasure) = transport_origin(parent(μ))
@inline from_origin(μ::PushforwardMeasure, x) = μ.f(from_origin(parent(μ), x))
@inline to_origin(μ::PushforwardMeasure, y) = μ.inv_f(to_origin(parent(μ), y))


function Base.rand(rng::AbstractRNG, ::Type{T}, ν::PushforwardMeasure) where {T}
    return from_origin(ν, rand(rng, T, transport_origin(ν)))
end

export pushfwd

"""
    pushfwd(f, μ, volcorr = WithVolCorr())

Return the [pushforward measure](https://en.wikipedia.org/wiki/Pushforward_measure)
from `μ` the [measurable function](https://en.wikipedia.org/wiki/Measurable_function) `f`.
"""
pushfwd(f, μ, volcorr = WithVolCorr()) = PushforwardMeasure(f, inverse(f), μ, volcorr)
