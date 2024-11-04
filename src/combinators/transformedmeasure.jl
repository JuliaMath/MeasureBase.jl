"""
    abstract type PushFwdStyle

Provides the behavior of a measure's [`rootmeasure`](@ref) under a
pushforward. Either [`AdaptRootMeasure()`](@ref) or
[`PushfwdRootMeasure()`](@ref)
"""
abstract type PushFwdStyle end
export PushFwdStyle

const TransformVolCorr = PushFwdStyle

"""
    AdaptRootMeasure()

Indicates that when applying a pushforward to a measure, it's
[`rootmeasure`](@ref) not not be pushed forward. Instead, the root measure
should be kept just "reshaped" to the new measurable space if necessary.

Density calculations for pushforward measures constructed with
`AdaptRootMeasure()` will take take the volume element of variate
transform (typically via the log-abs-det-Jacobian of the transform) into
account.
"""
struct AdaptRootMeasure <: TransformVolCorr end
export AdaptRootMeasure

const WithVolCorr = AdaptRootMeasure

"""
    PushfwdRootMeasure()

Indicates than when applying a pushforward to a measure, it's
[`rootmeasure`](@ref) should be pushed forward with the same function.

Density calculations for pushforward measures constructed with
`PushfwdRootMeasure()` will ignore the volume element of the variate
transform.
"""
struct PushfwdRootMeasure <: TransformVolCorr end
export PushfwdRootMeasure

const NoVolCorr = PushfwdRootMeasure

abstract type AbstractTransformedMeasure <: AbstractMeasure end

abstract type AbstractPushforward <: AbstractTransformedMeasure end

abstract type AbstractPullback <: AbstractTransformedMeasure end

function gettransform(::AbstractTransformedMeasure) end

function params(::AbstractTransformedMeasure) end

function paramnames(::AbstractTransformedMeasure) end

function parent(::AbstractTransformedMeasure) end

export PushforwardMeasure

"""
    struct PushforwardMeasure{F,I,M,S<:PushFwdStyle} <: AbstractPushforward
        f :: F
        finv :: I
        origin :: M
        style :: S
    end

    Users should not call `PushforwardMeasure` directly. Instead call or add
    methods to `pushfwd`.
"""
struct PushforwardMeasure{F,I,M,S<:PushFwdStyle} <: AbstractPushforward
    f::F
    finv::I
    origin::M
    style::S

    function PushforwardMeasure{F,I,M,S}(
        f::F,
        finv::I,
        origin::M,
        style::S,
    ) where {F,I,M,S<:PushFwdStyle}
        new{F,I,M,S}(f, finv, origin, style)
    end

    function PushforwardMeasure(f, finv, origin::M, style::S) where {M,S<:PushFwdStyle}
        new{Core.Typeof(f),Core.Typeof(finv),M,S}(f, finv, origin, style)
    end
end

const _NonBijectivePusfwdMeasure{M<:PushforwardMeasure,S<:PushFwdStyle} = Union{
    PushforwardMeasure{<:Any,<:NoInverse,M,S},
    PushforwardMeasure{<:NoInverse,<:Any,M,S},
    PushforwardMeasure{<:NoInverse,<:NoInverse,M,S},
}

gettransform(ν::PushforwardMeasure) = ν.f
parent(ν::PushforwardMeasure) = ν.origin

function Pretty.tile(ν::PushforwardMeasure)
    Pretty.list_layout(Pretty.tile.([ν.f, ν.origin]); prefix = :PushforwardMeasure)
end

# TODO: THIS IS ALMOST CERTAINLY WRONG 
# @inline function logdensity_rel(
#     ν::PushforwardMeasure{FF1,IF1,M1,<:AdaptRootMeasure},
#     β::PushforwardMeasure{FF2,IF2,M2,<:AdaptRootMeasure},
#     y,
# ) where {FF1,IF1,M1,FF2,IF2,M2}
#     x = β.inv_f(y)
#     f = ν.inv_f ∘ β.f
#     inv_f = β.inv_f ∘ ν.f
#     logdensity_rel(pushfwd(f, inv_f, ν.origin, AdaptRootMeasure()), β.origin, x)
# end

# TODO: Would profit from custom pullback:
function _combine_logd_with_ladj(logd_orig::Real, ladj::Real)
    logd_result = logd_orig + ladj
    R = typeof(logd_result)

    if isnan(logd_result) && isneginf(logd_orig) && isposinf(ladj)
        # Zero μ wins against infinite volume:
        R(-Inf)::R
    elseif isfinite(logd_orig) && isneginf(ladj)
        # Maybe  also for isneginf(logd_orig) && isfinite(ladj) ?
        # Return constant -Inf to prevent problems with ForwardDiff:
        #R(-Inf)
        near_neg_inf(R)::R # Avoids AdvancedHMC warnings
    else
        logd_result::R
    end
end

function logdensityof(
    @nospecialize(μ::_NonBijectivePusfwdMeasure{M,<:PushfwdRootMeasure}),
    @nospecialize(v::Any)
) where {M}
    throw(
        ArgumentError(
            "Can't calculate densities for non-bijective pushforward measure $(nameof(M))",
        ),
    )
end

function logdensityof(
    @nospecialize(μ::_NonBijectivePusfwdMeasure{M,<:AdaptRootMeasure}),
    @nospecialize(v::Any)
) where {M}
    throw(
        ArgumentError(
            "Can't calculate densities for non-bijective pushforward measure $(nameof(M))",
        ),
    )
end

for func in [:logdensityof, :logdensity_def]
    @eval function $func(ν::PushforwardMeasure{F,I,M,<:AdaptRootMeasure}, y) where {F,I,M}
        f_inv = unwrap(ν.finv)
        x, inv_ladj = with_logabsdet_jacobian(f_inv, y)
        logd_orig = $func(ν.origin, x)
        return _combine_logd_with_ladj(logd_orig, inv_ladj)
    end

    @eval function $func(ν::PushforwardMeasure{F,I,M,<:PushfwdRootMeasure}, y) where {F,I,M}
        f_inv = unwrap(ν.finv)
        x = f_inv(y)
        logd_orig = $func(ν.origin, x)
        return logd_orig
    end
end

insupport(m::PushforwardMeasure, x) = insupport(transport_origin(m), to_origin(m, x))

function testvalue(::Type{T}, ν::PushforwardMeasure) where {T}
    ν.f(testvalue(T, parent(ν)))
end

@inline function basemeasure(ν::PushforwardMeasure)
    pushfwd(ν.f, basemeasure(parent(ν)), PushfwdRootMeasure())
end

function rootmeasure(m::PushforwardMeasure{F,I,M,PushfwdRootMeasure}) where {F,I,M}
    pushfwd(m.f, rootmeasure(m.origin))
end
function rootmeasure(m::PushforwardMeasure{F,I,M,AdaptRootMeasure}) where {F,I,M}
    rootmeasure(m.origin)
end

_pushfwd_dof(::Type{MU}, ::Type, dof) where {MU} = NoDOF{MU}()
_pushfwd_dof(::Type{MU}, ::Type{<:Tuple{Any,Real}}, dof) where {MU} = dof

@inline getdof(ν::MU) where {MU<:PushforwardMeasure} = getdof(ν.origin)
@inline getdof(m::_NonBijectivePusfwdMeasure) = MeasureBase.NoDOF{typeof(m)}()

# Bypass `checked_arg`, would require potentially costly transformation:
@inline checked_arg(::PushforwardMeasure, x) = x

@inline transport_origin(ν::PushforwardMeasure) = ν.origin
@inline from_origin(ν::PushforwardMeasure, x) = ν.f(x)
@inline to_origin(ν::PushforwardMeasure, y) = ν.finv(y)

massof(m::PushforwardMeasure) = massof(transport_origin(m))

function Base.rand(rng::AbstractRNG, ::Type{T}, ν::PushforwardMeasure) where {T}
    return ν.f(rand(rng, T, ν.origin))
end

###############################################################################
# pushfwd

"""
    pushfwd(f, μ, style = AdaptRootMeasure())

Return the [pushforward
measure](https://en.wikipedia.org/wiki/Pushforward_measure) from `μ` the
[measurable function](https://en.wikipedia.org/wiki/Measurable_function) `f`.

To manually specify an inverse, call 
`pushfwd(InverseFunctions.setinverse(f, finv), μ, style)`.
"""
function pushfwd end
export pushfwd

@inline pushfwd(f, μ) = _pushfwd_impl(f, μ, AdaptRootMeasure())
@inline pushfwd(f, μ, style::AdaptRootMeasure) = _pushfwd_impl(f, μ, style)
@inline pushfwd(f, μ, style::PushfwdRootMeasure) = _pushfwd_impl(f, μ, style)

_pushfwd_impl(f, μ, style) = PushforwardMeasure(f, inverse(f), μ, style)

function _pushfwd_impl(
    f,
    μ::PushforwardMeasure{F,I,M,S},
    style::S,
) where {F,I,M,S<:PushFwdStyle}
    orig_μ = μ.origin
    new_f = fcomp(f, μ.f)
    new_f_inv = fcomp(μ.finv, inverse(f))
    PushforwardMeasure(new_f, new_f_inv, orig_μ, style)
end

_pushfwd_impl(::typeof(identity), μ, ::AdaptRootMeasure) = μ
_pushfwd_impl(::typeof(identity), μ, ::PushfwdRootMeasure) = μ

###############################################################################
# pullback

"""
    pullbck(f, μ, style = AdaptRootMeasure())

A _pullback_ is a dual concept to a _pushforward_. While a pushforward needs a
map _from_ the support of a measure, a pullback requires a map _into_ the
support of a measure. The log-density is then computed through function
composition, together with a volume correction as needed.

This can be useful, since the log-density of a `PushforwardMeasure` is computing
in terms of the inverse function; the "forward" function is not used at all. In
some cases, we may be focusing on log-density (and not, for example, sampling).

To manually specify an inverse, call 
`pullbck(InverseFunctions.setinverse(f, finv), μ, style)`.
"""
function pullbck end
export pullbck

@inline pullbck(f, μ) = _pullback_impl(f, μ, AdaptRootMeasure())
@inline pullbck(f, μ, style::AdaptRootMeasure) = _pullback_impl(f, μ, style)
@inline pullbck(f, μ, style::PushfwdRootMeasure) = _pullback_impl(f, μ, style)

function _pullback_impl(f, μ, style = AdaptRootMeasure())
    pushfwd(setinverse(inverse(f), f), μ, style)
end

@deprecate pullback(f, μ, style::PushFwdStyle = AdaptRootMeasure()) pullbck(f, μ, style)
