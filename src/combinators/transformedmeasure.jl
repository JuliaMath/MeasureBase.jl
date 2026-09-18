"""
    abstract type PushFwdStyle

Provides the behavior of a measure's [`rootmeasure`](@ref) under a
pushforward. Either [`AdaptRootMeasure()`](@ref) or
[`PushfwdRootMeasure()`](@ref)
"""
abstract type PushFwdStyle end
export PushFwdStyle

# Backward compatibility with user code, do not use in MeasureBase itself:
const TransformVolCorr = PushFwdStyle

"""
    AdaptRootMeasure()

Indicates that when applying a pushforward to a measure, its
[`rootmeasure`](@ref) need not be pushed forward. Instead, the root measure
should be kept just "reshaped" to the new measurable space if necessary.

Density calculations for pushforward measures constructed with
`AdaptRootMeasure()` will take take the volume element of variate
transform (typically via the log-abs-det-Jacobian of the transform) into
account.
"""
struct AdaptRootMeasure <: PushFwdStyle end
export AdaptRootMeasure

# Backward compatibility with user code, do not use in MeasureBase itself:
const WithVolCorr = AdaptRootMeasure

"""
    PushfwdRootMeasure()

Indicates that when applying a pushforward to a measure, its
[`rootmeasure`](@ref) should be pushed forward with the same function.

Density calculations for pushforward measures constructed with
`PushfwdRootMeasure()` will ignore the volume element of the variate
transform.
"""
struct PushfwdRootMeasure <: PushFwdStyle end
export PushfwdRootMeasure

# Backward compatibility with user code, do not use in MeasureBase itself:
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
struct PushforwardMeasure{F,I,M,S<:PushFwdStyle,VS} <: AbstractPushforward
    f::F
    finv::I
    origin::M
    style::S
    varsize::VS

    function PushforwardMeasure(f, finv, origin::M, style::S, varsize::VS) where {M,S<:PushFwdStyle,VS}
        new{Core.Typeof(f),Core.Typeof(finv),M,S,VS}(f, finv, origin, style, varsize)
    end
end

# The size of the variates of a pushforward follows from a test value of
# the origin, where the origin has variates of known size:
@inline function _pushfwd_varsize(f, μ)
    _pushfwd_varsize(f, μ, mspace_flatsize(μ))
end
@inline _pushfwd_varsize(f, μ, ::SizeLike) = _value_flatsize(f(testvalue(μ)))
@inline _pushfwd_varsize(f, μ, ::NoMSpaceElementSize) = NoMSpaceElementSize{typeof(μ)}()

@inline mspace_elsize(ν::PushforwardMeasure) = _value_or_unknown(ν.varsize, ν)
@inline mspace_flatsize(ν::PushforwardMeasure) = _value_or_unknown(ν.varsize, ν)
@inline _value_or_unknown(sz::SizeLike, ν) = sz
@inline _value_or_unknown(::NoMSpaceElementSize, ν) = NoMSpaceElementSize{typeof(ν)}()
@inline fixed_stream_size(::Type{<:PushforwardMeasure{<:Any,<:Any,<:Any,<:Any,VS}}) where {VS} = static(VS <: SizeLike)
@inline function mspace_ndims(::Type{MU}) where {VS,MU<:PushforwardMeasure{<:Any,<:Any,<:Any,<:Any,VS}}
    _ndims_of_size_type(VS, MU)
end
@inline _ndims_of_size_type(::Type{<:Tuple{Vararg{Any,N}}}, ::Type) where {N} = N
@inline _ndims_of_size_type(::Type{StaticArrays.Size{S}}, ::Type) where {S} = length(S)
@inline _ndims_of_size_type(::Type, ::Type{MU}) where {MU} = NoMSpaceElementSize{MU}()

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
function _combine_logd_with_ladj(logd_orig::Number, ladj::Number)
    logd_result = logd_orig + ladj
    R = typeof(logd_result)

    # Zero μ wins against infinite volume:
    zero_wins = isnan(logd_result) & isneginf(logd_orig) & isposinf(ladj)
    # Maybe also for isneginf(logd_orig) && isfinite(ladj) ?
    # Return near_neg_inf instead of constant -Inf to prevent problems
    # with ForwardDiff and to avoid AdvancedHMC warnings:
    fades_out = isfinite(logd_orig) & isneginf(ladj)

    ifelse(zero_wins, R(-Inf), ifelse(fades_out, near_neg_inf(R), logd_result))::R
end

function logdensityof_impl(
    @nospecialize(μ::_NonBijectivePusfwdMeasure{M,<:PushfwdRootMeasure}),
    @nospecialize(v::Any)
) where {M}
    throw(
        ArgumentError(
            "Can't calculate densities for non-bijective pushforward measure $(nameof(M))",
        ),
    )
end

function logdensityof_impl(
    @nospecialize(μ::_NonBijectivePusfwdMeasure{M,<:AdaptRootMeasure}),
    @nospecialize(v::Any)
) where {M}
    throw(
        ArgumentError(
            "Can't calculate densities for non-bijective pushforward measure $(nameof(M))",
        ),
    )
end

for (head, func) in [(:logdensityof_impl, :logdensityof), (:logdensity_def, :logdensity_def)]
    @eval function $head(ν::PushforwardMeasure{F,I,M,<:AdaptRootMeasure}, y) where {F,I,M}
        f_inv = unwrap(ν.finv)
        x, inv_ladj = with_logabsdet_jacobian(f_inv, y)
        logd_orig = $func(ν.origin, x)
        return _combine_logd_with_ladj(logd_orig, inv_ladj)
    end

    @eval function $head(ν::PushforwardMeasure{F,I,M,<:PushfwdRootMeasure}, y) where {F,I,M}
        f_inv = unwrap(ν.finv)
        x = f_inv(y)
        logd_orig = $func(ν.origin, x)
        return logd_orig
    end
end

# Checking insupport via the origin would require a potentially costly
# transformation of x:
insupport(m::PushforwardMeasure, x) = NoFastInsupport{typeof(m)}()

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

@inline fast_dof(ν::PushforwardMeasure) = fast_dof(ν.origin)
@inline fast_dof(m::_NonBijectivePusfwdMeasure) = MeasureBase.NoDOF{typeof(m)}()

# Bypass `checked_arg`, would require potentially costly transformation:
@inline checked_arg(::PushforwardMeasure, x) = x

# Pushforwards transport via their origin:
@inline transport_to_std(::Type{S}, ν::PushforwardMeasure, y) where {S<:StdMeasure} =
    transport_to_std(S, ν.origin, ν.finv(y))
@inline transport_from_std(::Type{S}, ν::PushforwardMeasure, z) where {S<:StdMeasure} =
    ν.f(transport_from_std(S, ν.origin, z))
@inline function transport_from_std_with_rest(::Type{S}, ν::PushforwardMeasure, z::AbstractVector) where {S<:StdMeasure}
    x, z_rest = transport_from_std_with_rest(S, ν.origin, z)
    return ν.f(x), z_rest
end

# Batched transport for pushforwards of measures with scalar variates, the
# functions apply elementwise then:
function batched_transport_to_std(::Type{S}, ν::PushforwardMeasure, Y::AbstractArray) where {S<:StdMeasure}
    _pushfwd_batched_to_std(S, ν, Y, mspace_flatsize(ν.origin))
end
@inline function _pushfwd_batched_to_std(::Type{S}, ν::PushforwardMeasure, Y::AbstractArray, ::Tuple{}) where {S}
    batched_transport_to_std(S, ν.origin, broadcast(ν.finv, Y))
end
@inline function _pushfwd_batched_to_std(::Type{S}, ν::PushforwardMeasure, Y::AbstractArray, ::Any) where {S}
    _batched_to_std(S, ν, Y, mspace_flatsize(ν))
end

function batched_transport_from_std(::Type{S}, ν::PushforwardMeasure, Z::AbstractArray) where {S<:StdMeasure}
    _pushfwd_batched_from_std(S, ν, Z, mspace_flatsize(ν.origin))
end
@inline function _pushfwd_batched_from_std(::Type{S}, ν::PushforwardMeasure, Z::AbstractArray, ::Tuple{}) where {S}
    broadcast(ν.f, batched_transport_from_std(S, ν.origin, Z))
end
@inline function _pushfwd_batched_from_std(::Type{S}, ν::PushforwardMeasure, Z::AbstractArray, ::Any) where {S}
    _batched_from_std(S, ν, Z, mspace_flatsize(ν))
end

massof(m::PushforwardMeasure) = massof(m.origin)

rand_impl(ctx::GenContext, ν::PushforwardMeasure) = ν.f(rand_impl(ctx, ν.origin))

# Batches of pushforwards apply the function to the variates of a batch of
# the origin, elementwise for scalar variates:
function batched_rand_impl(ctx::GenContext, ν::PushforwardMeasure, sz::Dims)
    _pushfwd_batched_rand(ctx, ν, sz, mspace_flatsize(ν.origin))
end
@inline function _pushfwd_batched_rand(ctx::GenContext, ν::PushforwardMeasure, sz::Dims, ::Tuple{})
    broadcast(ν.f, batched_rand_impl(ctx, ν.origin, sz))
end
@inline function _pushfwd_batched_rand(ctx::GenContext, ν::PushforwardMeasure, sz::Dims, ::Any)
    _batched_rand_pointwise(ctx, ν, sz)
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

@inline pushfwd(f) = Base.Fix1(pushfwd, f)
@inline pushfwd(f, μ) = _pushfwd_impl(f, μ, AdaptRootMeasure())
@inline pushfwd(f, μ, style::PushFwdStyle) = _pushfwd_impl(f, μ, style)

@inline pushfwd(::typeof(identity), μ) = μ
@inline pushfwd(::typeof(identity), μ, ::PushFwdStyle) = μ

_pushfwd_impl(f, μ, style) = PushforwardMeasure(f, inverse(f), μ, style, _pushfwd_varsize(f, μ))

function _pushfwd_impl(
    f,
    μ::PushforwardMeasure{F,I,M,S},
    style::S,
) where {F,I,M,S<:PushFwdStyle}
    orig_μ = μ.origin
    new_f = fcomp(f, μ.f)
    new_f_inv = fcomp(μ.finv, inverse(f))
    PushforwardMeasure(new_f, new_f_inv, orig_μ, style, _pushfwd_varsize(new_f, orig_μ))
end

# Simplifications for Dirac and WeightedMeasure origins are defined in
# smart-constructors.jl.

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

@inline pullbck(f) = Base.Fix1(pullbck, f)
@inline pullbck(f, μ) = _pullback_impl(f, μ, AdaptRootMeasure())
@inline pullbck(f, μ, style::PushFwdStyle) = _pullback_impl(f, μ, style)

function _pullback_impl(f, μ, style = AdaptRootMeasure())
    pushfwd(inverse(f), μ, style)
end

@deprecate pullback(f, μ, style::PushFwdStyle = AdaptRootMeasure()) pullbck(f, μ, style)
