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

Measures with counting root measures (e.g. point masses and discrete
distributions) are always pushed forward with
[`PushfwdRootMeasure()`](@ref): injective maps push counting measure
forward to counting measure on the image, without a volume element.
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

"""
    MeasureBase.gettransform(ν::AbstractTransformedMeasure)

Returns the function that `ν` transforms its origin measure with.
"""
function gettransform(::AbstractTransformedMeasure) end

paramnames(::AbstractTransformedMeasure) = ()

function parent(::AbstractTransformedMeasure) end

export PushforwardMeasure

"""
    struct PushforwardMeasure{F,I,M,S<:PushFwdStyle,VS} <: AbstractPushforward
        f :: F
        finv :: I
        origin :: M
        style :: S
        varsize :: VS
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

# Marks variate sizes of pushforwards that are not declared. They come from
# the origin for elementwise functions, else from the type `f` returns for
# variates of the origin's variate rank (by type inference, without
# evaluating `f`), or from `f` applied to a test value of the origin,
# computed on demand:
struct _UndeclaredVarSize end

@inline _declared_varsize(::Nothing) = _UndeclaredVarSize()
@inline _declared_varsize(sz::SizeLike) = canonical_size(sz)
@noinline _declared_varsize(sz) = throw(ArgumentError("The variate size of a pushforward must be a size tuple, got $sz"))

const _UndeclaredPushfwd{F,M} = PushforwardMeasure{F,<:Any,M,<:PushFwdStyle,_UndeclaredVarSize}

# Pushforwards by elementwise functions keep the variate size of their
# origin:
const _ElementwisePushfwd{M,S} = PushforwardMeasure{<:Base.BroadcastFunction,<:Base.BroadcastFunction,M,S}
const _ElementwiseUndeclaredPushfwd{M} =
    PushforwardMeasure{<:Base.BroadcastFunction,<:Base.BroadcastFunction,M,<:PushFwdStyle,_UndeclaredVarSize}

@inline mspace_flatsize(::Type{<:PushforwardMeasure{<:Any,<:Any,<:Any,<:Any,Tuple{}}}) = ()
@inline mspace_flatsize(::Type{<:PushforwardMeasure{<:Any,<:Any,<:Any,<:Any,StaticArrays.Size{S}}}) where {S} = StaticArrays.Size(S)
@inline mspace_flatsize(::Type{MU}) where {F,M,MU<:_UndeclaredPushfwd{F,M}} = _type_flatsize(_image_type(F, M), MU)
@inline mspace_flatsize(::Type{<:_ElementwiseUndeclaredPushfwd{M}}) where {M} = mspace_flatsize(M)

@inline function mspace_ndims(::Type{MU}) where {VS,MU<:PushforwardMeasure{<:Any,<:Any,<:Any,<:Any,VS}}
    _ndims_of_size_type(VS, MU)
end
@inline mspace_ndims(::Type{MU}) where {M,MU<:_ElementwisePushfwd{M}} = mspace_ndims(M)
@inline _ndims_of_size_type(::Type{<:Tuple{Vararg{Any,N}}}, ::Type) where {N} = N
@inline _ndims_of_size_type(::Type{StaticArrays.Size{S}}, ::Type) where {S} = length(S)
@inline _ndims_of_size_type(::Type{_UndeclaredVarSize}, ::Type{MU}) where {F,M,MU<:_UndeclaredPushfwd{F,M}} = _type_ndims(_image_type(F, M), MU)

@inline fixed_stream_size(::Type{<:PushforwardMeasure{<:Any,<:Any,<:Any,<:Any,VS}}) where {VS} = static(VS <: SizeLike)
@inline function fixed_stream_size(::Type{MU}) where {F,M,MU<:_UndeclaredPushfwd{F,M}}
    fixed_stream_size(M) & static(mspace_ndims(MU) isa Integer)
end

# The type `f` returns for variates of origins of type `M` with a fixed
# variate size, inferred for numbers resp. arrays of the variate rank (or
# static size) of `M`, `Union{}` where it isn't known:
@inline _image_type(::Type{F}, ::Type{M}) where {F,M} = _applied_type(F, _variate_standin(M, fixed_stream_size(M)))
@inline _applied_type(::Type{F}, ::Type{X}) where {F,X} = Core.Compiler.return_type(_apply_f, Tuple{F,X})
@inline _applied_type(::Type, ::Nothing) = Union{}
@inline _apply_f(f, x) = f(x)

@inline _variate_standin(::Type{M}, ::True) where {M} = _variate_standin(mspace_flatsize(M), mspace_ndims(M))
@inline _variate_standin(::Type, ::False) = nothing
@inline _variate_standin(::Tuple{}, ::Any) = Float64
@inline _variate_standin(::StaticArrays.Size{S}, ::Any) where {S} = StaticArrays.SArray{Tuple{S...},Float64,length(S),prod(S)}
@inline _variate_standin(::Union{SizeLike,NoMSpaceElementSize}, n) = _rank_standin(_static_ndims_of(n))
@inline _rank_standin(::StaticInteger{0}) = Float64
@inline _rank_standin(::StaticInteger{N}) where {N} = Array{Float64,N}
@inline _rank_standin(::NoMSpaceElementSize) = nothing

@inline function _type_flatsize(::Type{T}, ::Type{MU}) where {T,MU}
    if T !== Union{} && T <: Number
        ()
    elseif T !== Union{} && T <: StaticArray && T isa DataType && eltype(T) <: Number
        StaticArrays.Size(T)
    else
        NoMSpaceElementSize{MU}()
    end
end
@inline function _type_ndims(::Type{T}, ::Type{MU}) where {T,MU}
    if T !== Union{} && T <: Number
        0
    elseif T !== Union{} && T <: AbstractArray && T isa DataType && eltype(T) <: Number
        ndims(T)
    else
        NoMSpaceElementSize{MU}()
    end
end

# The image of a test value of the origin:
_testvalue_image(ν::PushforwardMeasure) = ν.f(testvalue(ν.origin))

# Flat variate sizes known without evaluating `f`:
@inline _known_flatsize(::PushforwardMeasure, sz::SizeLike) = sz
@inline _known_flatsize(ν::_ElementwisePushfwd, ::_UndeclaredVarSize) = mspace_flatsize(ν.origin)
@inline _known_flatsize(ν::PushforwardMeasure, ::_UndeclaredVarSize) = mspace_flatsize(typeof(ν))

# Other sizes are learned from the image of a test value, for origins with
# variates of fixed size (recomputed on every call):
@inline mspace_flatsize(ν::PushforwardMeasure) = _value_or_unknown(_learn_flatsize(ν, _known_flatsize(ν, ν.varsize)), ν)
@inline mspace_elsize(ν::PushforwardMeasure) = mspace_flatsize(ν)
@inline _learn_flatsize(::PushforwardMeasure, sz::SizeLike) = sz
@inline _learn_flatsize(ν::PushforwardMeasure, ::NoMSpaceElementSize) = _learn_flatsize(ν, fixed_stream_size(ν.origin))
@inline _learn_flatsize(ν::PushforwardMeasure, ::True) = _pushfwd_varsize(ν)
@inline _learn_flatsize(ν::PushforwardMeasure, ::False) = NoMSpaceElementSize{typeof(ν)}()
@inline _value_or_unknown(sz::SizeLike, ν) = sz
@inline _value_or_unknown(::NoMSpaceElementSize, ν) = NoMSpaceElementSize{typeof(ν)}()

# Learning the size never takes part in automatic differentiation:
_pushfwd_varsize(ν::PushforwardMeasure) = _adignore_call(() -> _value_flatsize(_checked_testvalue_image(ν)))

function _checked_testvalue_image(ν::PushforwardMeasure)
    try
        _testvalue_image(ν)
    catch err
        _throw_unlearnable_varsize(ν, err)
    end
end
@noinline function _throw_unlearnable_varsize(ν, err)
    throw(ArgumentError("The variate size of a pushforward measure of type $(nameof(typeof(ν))) is learned from the image of a test value of its origin, which failed ($(sprint(showerror, err))). Declare the size via `pushfwd(f, μ; varsize)`."))
end

# Point densities need no variate size, the origin checks the variate
# after the inverse:
@inline function _point_ld(f::F, ν::PushforwardMeasure, x::AbstractArray{<:Number}) where {F}
    _pushfwd_point_ld(f, ν, x, _known_flatsize(ν, ν.varsize))
end
@inline _pushfwd_point_ld(f::F, ν, x, ::SizeLike) where {F} = _point_ld_byrank(f, ν, x, _static_ndims(ν))
@inline _pushfwd_point_ld(f::F, ν, x, ::NoMSpaceElementSize) where {F} = f(ν, x)

const _NonBijectivePusfwdMeasure{M,S<:PushFwdStyle} = Union{
    PushforwardMeasure{<:Any,<:NoInverse,M,S},
    PushforwardMeasure{<:NoInverse,<:Any,M,S},
    PushforwardMeasure{<:NoInverse,<:NoInverse,M,S},
}

gettransform(ν::PushforwardMeasure) = ν.f
parent(ν::PushforwardMeasure) = ν.origin

function Pretty.tile(ν::PushforwardMeasure)
    Pretty.list_layout(Pretty.tile.([ν.f, ν.origin]); prefix = :PushforwardMeasure)
end

# TODO: Would profit from custom pullback:
function _combine_logd_with_ladj(logd_orig::Number, ladj::Number)
    logd_result = logd_orig + ladj
    # Zero μ wins against infinite volume:
    zero_wins = isnan(logd_result) & isneginf(logd_orig) & isposinf(ladj)
    ifelse(zero_wins, oftype(logd_result, -Inf), logd_result)
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

# Pushforwards that adapt the root measure have the root measure of their
# origin as their base measure, their density relative to it is the full
# density of the origin with the volume element:
@inline function logdensityof_impl(ν::PushforwardMeasure{F,I,M,<:AdaptRootMeasure}, y) where {F,I,M}
    x, inv_ladj = with_logabsdet_jacobian(unwrap(ν.finv), y)
    return _combine_logd_with_ladj(logdensityof(ν.origin, x), inv_ladj)
end
@inline logdensity_def(ν::PushforwardMeasure{F,I,M,<:AdaptRootMeasure}, y) where {F,I,M} = logdensityof_impl(ν, y)

# Pushforwards that push the root measure forward have the pushforward of
# the base measure of their origin as their base measure:
for (head, func) in [(:logdensityof_impl, :logdensityof), (:logdensity_def, :logdensity_def)]
    @eval function $head(ν::PushforwardMeasure{F,I,M,<:PushfwdRootMeasure}, y) where {F,I,M}
        return $func(ν.origin, unwrap(ν.finv)(y))
    end
end

# Pushforwards by elementwise functions evaluate densities over flat
# batches, the log-abs-det-Jacobian terms sum over the variate dimensions
# of the origin:
for bhead in (:batched_logdensityof_impl, :batched_logdensity_def)
    @eval function $bhead(ν::_ElementwisePushfwd{M,<:AdaptRootMeasure}, Y) where {M}
        _elementwise_pushfwd_ld(logdensityof_impl, ν, Y, _static_ndims(ν.origin))
    end
end
for (bhead, head) in [(:batched_logdensityof_impl, :logdensityof_impl), (:batched_logdensity_def, :logdensity_def)]
    @eval function $bhead(ν::_ElementwisePushfwd{M,<:PushfwdRootMeasure}, Y) where {M}
        _batched_kernel($head, ν.origin, broadcast(ν.finv.f, Y))
    end
end

function _elementwise_pushfwd_ld(f::F, ν::PushforwardMeasure, Y, k::StaticInteger) where {F}
    f_inv = ν.finv.f
    ℓ = _batched_kernel(f, ν.origin, broadcast(f_inv, Y))
    ladj = sum_leading_dims(broadcast(_LadjOf(f_inv), Y), k)
    return _lazy_combine_ladj(ℓ, ladj)
end
function _elementwise_pushfwd_ld(f::F, ν::PushforwardMeasure, Y, ::NoMSpaceElementSize) where {F}
    _default_batched_kernel(f, ν, Y, _static_ndims(ν))
end

struct _LadjOf{F} <: Function
    f::F
end
@inline (k::_LadjOf)(y) = last(with_logabsdet_jacobian(k.f, y))

@inline _lazy_combine_ladj(ℓ::Number, ladj::Number) = _combine_logd_with_ladj(ℓ, ladj)
@inline _lazy_combine_ladj(ℓ, ladj) = Broadcast.instantiate(Broadcast.broadcasted(_combine_logd_with_ladj, ℓ, ladj))

# Checking the support via the origin would require a potentially costly
# transformation of x:
insupportof(m::PushforwardMeasure, x) = NoFastInsupport{typeof(m)}()

function testvalue(::Type{T}, ν::PushforwardMeasure) where {T}
    ν.f(testvalue(T, parent(ν)))
end

@inline function basemeasure(ν::PushforwardMeasure{F,I,M,PushfwdRootMeasure}) where {F,I,M}
    _pushfwd_base(ν, basemeasure(ν.origin))
end
# Same function, same variate size:
@inline _pushfwd_base(ν::PushforwardMeasure, β) = PushforwardMeasure(ν.f, ν.finv, β, ν.style, ν.varsize)

# The root measure of the origin, reshaped where the pushforward changes
# the variate rank:
@inline function basemeasure(ν::PushforwardMeasure{F,I,M,AdaptRootMeasure}) where {F,I,M}
    _reshape_root(rootmeasure(ν.origin), ν)
end
@inline _reshape_root(root, ::PushforwardMeasure) = root

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

# Batches of pushforwards apply the functions to flat batches of the
# origin, elementwise for `Base.BroadcastFunction`s and variate by variate
# (in a host loop) otherwise. The AffineMaps extension adds affine maps.
function batched_transport_to_std(::Type{S}, ν::PushforwardMeasure, Y) where {S<:StdMeasure}
    batched_transport_to_std(S, ν.origin, _apply_batched(ν.finv, Y, _static_ndims(ν)))
end
function batched_transport_from_std(::Type{S}, ν::PushforwardMeasure, Z::AbstractArray) where {S<:StdMeasure}
    _apply_batched(ν.f, batched_transport_from_std(S, ν.origin, Z), _static_ndims(ν.origin))
end

# Apply `f` to a flat batch of variates of rank `k`:
@inline _apply_batched(f, X, k) = _apply_generic(unwrap(f), X, k)
@inline _apply_generic(f, X, k) = _apply_by_rank(f, X, k)
@inline _apply_generic(f::Base.BroadcastFunction, X, k) = broadcast(f.f, X)
@inline _apply_by_rank(f, X, ::StaticInteger{0}) = broadcast(f, X)
@inline _apply_by_rank(f, X::AbstractArray, ::StaticInteger{0}) = broadcast(f, X)
@inline _apply_by_rank(f, X::AbstractArray, ::StaticInteger{K}) where {K} = _apply_to_slices(f, X, Val(K))
@inline _apply_to_slices(f, X::AbstractArray{<:Any,K}, ::Val{K}) where {K} = f(X)
@inline _apply_to_slices(f, X::AbstractArray, ::Val{K}) where {K} = stacked(map(f, sliced(X, Val(K))))
@noinline function _apply_by_rank(f, X, ::NoMSpaceElementSize)
    throw(ArgumentError("Applying functions of type $(nameof(typeof(f))) to batches of variates requires MeasureBase.mspace_ndims to be declared for the measure"))
end

massof(m::PushforwardMeasure) = massof(m.origin)
@inline _unweighted(ν::PushforwardMeasure) = PushforwardMeasure(ν.f, ν.finv, _unweighted(ν.origin), ν.style, ν.varsize)

rand_impl(ctx::GenContext, ν::PushforwardMeasure) = _pushfwd_draw(ctx, ν.origin, ν.f(rand_impl(ctx, ν.origin)))

# Transformed draws have the precision of the context, transformed atoms
# stay exact:
@inline _pushfwd_draw(ctx::GenContext, ::Any, y) = convert_realtype(get_precision(ctx), y)
@inline _pushfwd_draw(::GenContext, ::Dirac, y) = y

# Batches of pushforwards apply the function to the variates of a batch of
# the origin:
function batched_rand_impl(ctx::GenContext, ν::PushforwardMeasure, sz::SizeLike)
    _pushfwd_batched_rand(ctx, ν, sz, _static_ndims(ν.origin))
end
@inline function _pushfwd_batched_rand(ctx::GenContext, ν::PushforwardMeasure, sz::SizeLike, k::StaticInteger)
    _pushfwd_draw(ctx, ν.origin, _apply_batched(ν.f, batched_rand_impl(ctx, ν.origin, sz), k))
end
@inline function _pushfwd_batched_rand(ctx::GenContext, ν::PushforwardMeasure, sz::SizeLike, ::NoMSpaceElementSize)
    _batched_rand_pointwise(ctx, ν, sz)
end

###############################################################################
# pushfwd

"""
    pushfwd(f, μ, style = AdaptRootMeasure(); varsize)

Return the [pushforward
measure](https://en.wikipedia.org/wiki/Pushforward_measure) of `μ` under the
[measurable function](https://en.wikipedia.org/wiki/Measurable_function) `f`.

To manually specify an inverse, call
`pushfwd(InverseFunctions.setinverse(f, finv), μ, style)`.

Pushforwards are only defined on the image of `f`. Kernels of measures
never throw for values in their variate space, but the inverse of `f` may
throw outside of the image: `logdensityof(pushfwd(exp, StdNormal()), -1.0)`
throws a `DomainError` from `log`.

`varsize` declares the size of the variates of the pushforward (`()` for
scalars). Without it, the size is that of `μ` for
`Base.BroadcastFunction`s, follows from the type `f` returns for variates
of the variate rank (or static size) of `μ`, or else is learned from `f`
applied to a test value of `μ` (see [`testvalue`](@ref)) when an
operation needs it (batched kernels, streams). `f` is never evaluated at
construction, so `μ` needs no random variates unless the size must be
learned.
"""
function pushfwd end
export pushfwd

@inline pushfwd(f) = Base.Fix1(pushfwd, f)
@inline pushfwd(f, μ; varsize = nothing) = _pushfwd_impl(f, μ, AdaptRootMeasure(), _declared_varsize(varsize))
@inline pushfwd(f, μ, style::PushFwdStyle; varsize = nothing) = _pushfwd_impl(f, μ, style, _declared_varsize(varsize))

@inline pushfwd(::typeof(identity), μ; varsize = nothing) = μ
@inline pushfwd(::typeof(identity), μ, ::PushFwdStyle; varsize = nothing) = μ

_pushfwd_impl(f, μ, style, varsize) = _pushfwd_with_style(f, μ, _pushfwd_style(style, μ), varsize)

_pushfwd_with_style(f, μ, style, varsize) = PushforwardMeasure(f, inverse(f), μ, style, varsize)

function _pushfwd_with_style(
    f,
    μ::PushforwardMeasure{F,I,M,S},
    style::S,
    varsize,
) where {F,I,M,S<:PushFwdStyle}
    PushforwardMeasure(fcomp(f, μ.f), fcomp(μ.finv, inverse(f)), μ.origin, style, varsize)
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

function Adapt.adapt_structure(to, ν::PushforwardMeasure)
    PushforwardMeasure(Adapt.adapt(to, ν.f), Adapt.adapt(to, ν.finv), Adapt.adapt(to, ν.origin), ν.style, ν.varsize)
end
