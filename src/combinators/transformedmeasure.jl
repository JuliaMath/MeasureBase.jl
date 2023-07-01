# TODO: Compare with ChangesOfVariables.jl

using InverseFunctions: FunctionWithInverse

abstract type AbstractTransformedMeasure <: AbstractMeasure end

abstract type AbstractPushforward <: AbstractTransformedMeasure end

abstract type AbstractPullback <: AbstractTransformedMeasure end

function gettransform(::AbstractTransformedMeasure) end

function params(::AbstractTransformedMeasure) end

function paramnames(::AbstractTransformedMeasure) end

function parent(::AbstractTransformedMeasure) end

export PushforwardMeasure

"""
    struct PushforwardMeasure{F,I,M,VC<:TransformVolCorr} <: AbstractPushforward
        f :: F
        finv :: I
        origin :: M
        volcorr :: VC
    end

    Users should not call `PushforwardMeasure` directly. Instead call or add
    methods to `pushfwd`.
"""
struct PushforwardMeasure{F,I,M,VC<:TransformVolCorr} <: AbstractPushforward
    f::F
    finv::I
    origin::M
    volcorr::VC
end

gettransform(ν::PushforwardMeasure) = ν.f
parent(ν::PushforwardMeasure) = ν.origin

function Pretty.tile(ν::PushforwardMeasure)
    Pretty.list_layout(Pretty.tile.([ν.f, ν.origin]); prefix = :PushforwardMeasure)
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

@inline function logdensity_def(ν::PushforwardMeasure{F,I,M,<:WithVolCorr}, y) where {F,I,M}
    f = ν.f
    finv = ν.finv
    x_orig, inv_ladj = with_logabsdet_jacobian(unwrap(finv), y)
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

@inline function logdensity_def(ν::PushforwardMeasure{F,I,M,<:NoVolCorr}, y) where {F,I,M}
    x = ν.finv(y)
    return logdensity_def(ν.origin, x)
end

insupport(ν::PushforwardMeasure, y) = insupport(ν.origin, ν.finv(y))

function testvalue(::Type{T}, ν::PushforwardMeasure) where {T}
    ν.f(testvalue(T, parent(ν)))
end

@inline function basemeasure(ν::PushforwardMeasure)
    pushfwd(ν.f, basemeasure(parent(ν)), NoVolCorr())
end


const _NonBijectivePushforward = Union{PushforwardMeasure{<:Any,<:NoInverse},PushforwardMeasure{<:NoInverse,<:Any},PushforwardMeasure{<:NoInverse,<:NoInverse}}

@inline getdof(ν::PushforwardMeasure) = _pushfwd_dof(ν)
_pushfwd_dof(ν::PushforwardMeasure) = getdof(ν.origin)
_pushfwd_dof(ν::_NonBijectivePushforward) = NoDOF{typeof(ν)}()

@inline fast_dof(ν::PushforwardMeasure) = _pushfwd_fastdof(ν)
_pushfwd_fastdof(ν::PushforwardMeasure) = fast_dof(ν.origin)
_pushfwd_fastdof(ν::_NonBijectivePushforward) = NoDOF{typeof(ν)}()


# Bypass `checked_arg`, would require potentially costly transformation:
@inline checked_arg(::PushforwardMeasure, x) = x

@inline transport_origin(ν::PushforwardMeasure) = ν.origin
@inline from_origin(ν::PushforwardMeasure, x) = ν.f(x)
@inline to_origin(ν::PushforwardMeasure, y) = ν.finv(y)

function Base.rand(rng::AbstractRNG, ::Type{T}, ν::PushforwardMeasure) where {T}
    return ν.f(rand(rng, T, parent(ν)))
end

###############################################################################
# pushfwd

export pushfwd

"""
    pushfwd(f, μ, volcorr = WithVolCorr())

Return the [pushforward
measure](https://en.wikipedia.org/wiki/Pushforward_measure) from `μ` the
[measurable function](https://en.wikipedia.org/wiki/Measurable_function) `f`.

To manually specify an inverse, call 
`pushfwd(InverseFunctions.setinverse(f, finv), μ, volcorr)`.
"""
function pushfwd(f, μ, volcorr::TransformVolCorr = WithVolCorr())
    PushforwardMeasure(f, inverse(f), μ, volcorr)
end

function pushfwd(f, μ::PushforwardMeasure, volcorr::TransformVolCorr = WithVolCorr())
    _pushfwd_of_pushfwd(f, μ, μ.volcorr, volcorr)
end

# Either both WithVolCorr or both NoVolCorr, so we can merge them
function _pushfwd_of_pushfwd(f, μ::PushforwardMeasure, ::V, v::V) where {V}
    pushfwd(fchain((μ.f, f)), μ.origin, v)
end

function _pushfwd_of_pushfwd(f, μ::PushforwardMeasure, _, v)
    PushforwardMeasure(f, inverse(f), μ, v)
end

###############################################################################
# pullback

"""
    pullbck(f, μ, volcorr = WithVolCorr())

A _pullback_ is a dual concept to a _pushforward_. While a pushforward needs a
map _from_ the support of a measure, a pullback requires a map _into_ the
support of a measure. The log-density is then computed through function
composition, together with a volume correction as needed.

This can be useful, since the log-density of a `PushforwardMeasure` is computing
in terms of the inverse function; the "forward" function is not used at all. In
some cases, we may be focusing on log-density (and not, for example, sampling).

To manually specify an inverse, call 
`pullbck(InverseFunctions.setinverse(f, finv), μ, volcorr)`.
"""
function pullbck(f, μ, volcorr::TransformVolCorr = WithVolCorr())
    PushforwardMeasure(inverse(f), f, μ, volcorr)
end
export pullbck

@deprecate pullback(f, μ, volcorr::TransformVolCorr = WithVolCorr()) pullbck(f, μ, volcorr)
