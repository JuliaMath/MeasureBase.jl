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
    struct PushforwardMeasure{F,MU,VC<:TransformVolCorr} <: AbstractPushforward
        f :: F
        origin :: MU
        volcorr :: VC
    end
"""
struct PushforwardMeasure{F,M,VC<:TransformVolCorr} <: AbstractPushforward
    f::F
    origin::M
    volcorr::VC
end

gettransform(ν::PushforwardMeasure) = ν.f
parent(ν::PushforwardMeasure) = ν.origin

function transport_def(ν::PushforwardMeasure{F,M}, μ::M, x) where {F,M}
    if μ == parent(ν)
        return ν.f(x)
    else
        invoke(transport_def, Tuple{PushforwardMeasure,Any,Any}, ν, μ, x)
    end
end

function transport_def(μ::M, ν::PushforwardMeasure{F,M}, y) where {F,M}
    if μ == parent(ν)
        return inverse(ν.f)(y)
    else
        invoke(transport_def, Tuple{Any,PushforwardMeasure,Any}, μ, ν, y)
    end
end

function Pretty.tile(ν::PushforwardMeasure)
    Pretty.list_layout(Pretty.tile.([ν.f, ν.origin]); prefix = :PushforwardMeasure)
end

# TODO: Reduce code duplication
@inline function logdensityof(
    ν::PushforwardMeasure{F,M,<:WithVolCorr},
    y,
) where {F,M}
    f = ν.f 
    finv = inverse(f)
    x_orig, inv_ladj = with_logabsdet_jacobian(finv, y)
    μ = ν.origin
    logd_orig = unsafe_logdensityof(ν.origin, x_orig)
    logd = float(logd_orig + inv_ladj)
    neginf = oftype(logd, -Inf)
    insupport(μ, x_orig) || return oftype(logd, -Inf)
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

# TODO: THIS IS ALMOST CERTAINLY WRONG 
# @inline function logdensity_rel(
#     ν::PushforwardMeasure{FF1,IF1,M1,<:WithVolCorr},
#     β::PushforwardMeasure{FF2,IF2,M2,<:WithVolCorr},
#     y,
# ) where {FF1,IF1,M1,FF2,IF2,M2}
#     x = β.inv_f(y)
#     f = ν.inv_f ∘ β.f
#     inv_f = β.inv_f ∘ ν.f
#     logdensity_rel(pushfwd(f, inv_f, ν.origin, WithVolCorr()), β.origin, x)
# end

@inline function logdensity_def(
    ν::PushforwardMeasure{F,M,<:WithVolCorr},
    y,
) where {F,M}
    f = ν.f
    finv = inverse(f)
    x_orig, inv_ladj = with_logabsdet_jacobian(finv, y)
    logd_orig = unsafe_logdensityof(ν.origin, x_orig)
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
    ν::PushforwardMeasure{F,M,<:NoVolCorr},
    y,
) where {F,M}
    x_orig = to_origin(ν, y)
    return unsafe_logdensityof(ν.origin, x_orig)
end

insupport(ν::PushforwardMeasure, y) = insupport(transport_origin(ν), to_origin(ν, y))

function testvalue(::Type{T}, ν::PushforwardMeasure) where {T}
    from_origin(ν, testvalue(T, transport_origin(ν)))
end

@inline function basemeasure(ν::PushforwardMeasure)
    PushforwardMeasure(ν.f, basemeasure(transport_origin(ν)), NoVolCorr())
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

@inline transport_origin(ν::PushforwardMeasure) = ν.origin
@inline from_origin(ν::PushforwardMeasure, x) = ν.f(x)
@inline to_origin(ν::PushforwardMeasure, y) = inverse(ν.f)(y)

function Base.rand(rng::AbstractRNG, ::Type{T}, ν::PushforwardMeasure) where {T}
    return ν.f(rand(rng, T, parent(ν)))
end

export pushfwd

"""
    pushfwd(f, [f_inverse,] μ, volcorr = WithVolCorr())

Return the [pushforward
measure](https://en.wikipedia.org/wiki/Pushforward_measure) from `μ` the
[measurable function](https://en.wikipedia.org/wiki/Measurable_function) `f`.

If `f_inverse` is specified, it must be a valid inverse of the function given by
restricting `f` to the support of `μ`.
"""
function pushfwd(f, μ::AbstractMeasure, volcorr::TransformVolCorr = WithVolCorr())
    PushforwardMeasure(f, μ, volcorr)
end

function pushfwd(f, finv, μ::AbstractMeasure, volcorr::TransformVolCorr = WithVolCorr())
    pushfwd(setinverse(f, finv), μ, volcorr)
end

"""
    pullback(f, [f_inverse,] μ, volcorr = WithVolCorr())

A _pullback_ is a dual concept to a _pushforward_. While a pushforward needs a
map _from_ the support of a measure, a pullback requires a map _into_ the
support of a measure. The log-density is then computed through function
composition, together with a volume correction as needed.

This can be useful, since the log-density of a `PushforwardMeasure` is computing
in terms of the inverse function; the "forward" function is not used at all. In
some cases, we may be focusing on log-density (and not, for example, sampling).
"""
function pullback(f, μ::AbstractMeasure, volcorr::TransformVolCorr = WithVolCorr())
    pushfwd(inverse(f), μ, volcorr)
end

function pullback(f, finv, μ::AbstractMeasure, volcorr::TransformVolCorr = WithVolCorr())
    pushfwd(setinverse(finv, f), μ, volcorr)
end

# TODO: implement pushfwd-of-a-pushforward
function pushfwd(f, μ::PushforwardMeasure, volcorr::TransformVolCorr = WithVolCorr())
    _pushfwd(f, μ, μ.volcorr, volcorr)
end

# Either both WithVolCorr or both NoVolCorr, so we can merge them
function _pushfwd(f, μ, ::V, v::V) where {V}
    pushfwd(f ∘ μ.f, μ.origin, v)
end

function _pushfwd(f, μ, _, v)
    PushforwardMeasure(f, μ, v)
end
