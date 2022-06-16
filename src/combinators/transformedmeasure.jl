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
    struct PushforwardMeasure{FF,IF,M} <: AbstractPushforward
        f :: FF
        inv_f :: IF
        origin :: M
    end
"""
struct PushforwardMeasure{FF,IF,M} <: AbstractPushforward
    f::FF
    inv_f::IF
    origin::M
end

gettransform(μ::PushforwardMeasure) = μ.f
parent(μ::PushforwardMeasure) = μ.origin


function Pretty.tile(μ::PushforwardMeasure)
    Pretty.list_layout(Pretty.tile.([μ.f, μ.inv_f, μ.origin]); prefix = :PushforwardMeasure)
end

@inline function logdensity_def(μ::PushforwardMeasure, x)
    x_orig, inv_ladj = with_logabsdet_jacobian(μ.inv_f, x)
    logd_orig = logdensityof(μ.origin, x_orig)

    logd = logd_orig + inv_ladj
    R = typeof(logd)
    if isnan(logd) && logd_orig == -Inf && inv_ladj == +Inf
        # Zero density wins against infinite volume:
        R(-Inf)
    elseif isfinite(logd_orig) && (inv_ladj == -Inf)
        # Maybe  also for (logd_orig == -Inf) && isfinite(inv_ladj) ?
        # Return constant -Inf to prevent problems with ForwardDiff:
        R(-Inf)
    else
        logd
    end
end


insupport(μ::PushforwardMeasure, x) = insupport(to_origin(μ, x))

testvalue(μ::PushforwardMeasure) = from_origin(μ, testvalue(vartransform_origin(μ)))

@inline function basemeasure(μ::PushforwardMeasure)
    PushforwardMeasure(μ.f, μ.inv_f, basemeasure(vartransform_origin(μ)))
end

@inline getdof(::MU) where {MU<:PushforwardMeasure} = NoDOF{MU}()

@inline checked_var(::MU, ::Any) where {MU<:PushforwardMeasure} = NoVarCheck{MU}()

@inline vartransform_origin(μ::PushforwardMeasure) = μ.origin
@inline to_origin(μ::PushforwardMeasure, y) = μ.inv_f(y)
@inline from_origin(μ::PushforwardMeasure, x) = μ.f(x)

function Base.rand(rng::AbstractRNG, ::Type{T}, μ::PushforwardMeasure) where T
    return from_origin(μ, rand(rng, T, vartransform_origin(μ)))
end


export pushfwd

"""
    pushfwd(f, μ)

Return the [pushforward measure](https://en.wikipedia.org/wiki/Pushforward_measure)
from `μ` the [measurable function](https://en.wikipedia.org/wiki/Measurable_function) `f`.
"""
pushfwd(f, μ) = PushforwardMeasure(f, inverse(f), μ)
