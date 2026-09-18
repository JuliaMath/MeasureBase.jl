import Base
import StaticThings

export PowerMeasure

"""
    struct PowerMeasure{M,...} <: AbstractProductMeasure

A power measure is a product of a measure with itself. The number of elements in
the product determines the dimensionality of the resulting support.

Note that power measures are only well-defined for integer powers.

The nth power of a measure μ can be written μ^n.

See also [`pwr_base`](@ref), [`pwr_axes`](@ref) and [`pwr_size`](@ref).
"""
struct PowerMeasure{M,A} <: AbstractProductMeasure
    parent::M
    axes::A
end

StaticThings.maybestatic_length(μ::PowerMeasure) = size2length(maybestatic_size(μ))
StaticThings.maybestatic_size(μ::PowerMeasure) = axes2size(μ.axes)

"""
    MeasureBase.pwr_base(μ::PowerMeasure)

Returns `ν` for `μ = ν^axs`
"""
@inline pwr_base(μ::PowerMeasure) = μ.parent

"""
    MeasureBase.pwr_axes(μ::PowerMeasure)

Returns `axs` for `μ = ν^axs`, `axs` being a tuple of integer ranges.
"""
@inline pwr_axes(μ::PowerMeasure) = μ.axes

"""
    MeasureBase.pwr_size(μ::PowerMeasure)

Returns `sz` for `μ = ν^sz`, `sz` being a tuple of integers.
"""
@inline pwr_size(μ::PowerMeasure) = axes2size(μ.axes)

function Pretty.tile(μ::PowerMeasure)
    sz = length.(μ.axes)
    arg1 = Pretty.tile(μ.parent)
    arg2 = Pretty.tile(length(sz) == 1 ? only(sz) : sz)
    return Pretty.pair_layout(arg1, arg2; sep = " ^ ")
end

function _cartidxs(axs::Tuple{Vararg{AbstractUnitRange,N}}) where {N}
    CartesianIndices(map(asnonstatic, axs))
end

# Variates of powers are generated as one flat batch of variates of the
# base measure, with the power's size as additional batch dimensions. Base
# measures without fixed variate sizes generate their variates one by one.

rand_impl(ctx::GenContext, μ::PowerMeasure) = _pwr_rand(ctx, μ, fixed_stream_size(pwr_base(μ)))
_pwr_rand(ctx::GenContext, μ::PowerMeasure, ::True) = _pwr_variate(μ, batched_rand_impl(ctx, μ, ()))
function _pwr_rand(ctx::GenContext, μ::PowerMeasure, ::False)
    ν = pwr_base(μ)
    map(_ -> rand_impl(ctx, ν), _cartidxs(pwr_axes(μ)))
end

function batched_rand_impl(ctx::GenContext, μ::PowerMeasure, sz::Dims)
    _pwr_batched_rand(ctx, μ, sz, fixed_stream_size(pwr_base(μ)))
end
function _pwr_batched_rand(ctx::GenContext, μ::PowerMeasure, sz::Dims, ::True)
    batched_rand_impl(ctx, pwr_base(μ), (_dynamic_dims(pwr_size(μ))..., sz...))
end
_pwr_batched_rand(ctx::GenContext, μ::PowerMeasure, sz::Dims, ::False) = _batched_rand_pointwise(ctx, μ, sz)

marginals(d::PowerMeasure) = maybestatic_fill(d.parent, d.axes)

@inline mspace_elsize(μ::PowerMeasure) = pwr_size(μ)
@inline mspace_flatsize(μ::PowerMeasure) = _cat_sizes(mspace_flatsize(pwr_base(μ)), pwr_size(μ))
@inline function mspace_flatsize(::Type{<:PowerMeasure{M,A}}) where {M,A<:Tuple{Vararg{StaticOneToLike}}}
    _cat_sizes(mspace_flatsize(M), _static_axes_size(A))
end
@generated function _static_axes_size(::Type{A}) where {A<:Tuple{Vararg{StaticOneToLike}}}
    :(StaticArrays.Size($(map(T -> T.parameters[1], A.parameters)...)))
end

function Base.:^(μ::AbstractMeasure, dims::Tuple{Vararg{AbstractArray,N}}) where {N}
    powermeasure(μ, dims)
end

Base.:^(μ::AbstractMeasure, dims::Tuple) = powermeasure(μ, maybestatic_oneto.(dims))
Base.:^(μ::AbstractMeasure, n) = powermeasure(μ, (n,))

# Base.show(io::IO, d::PowerMeasure) = print(io, d.parent, " ^ ", size(d.xs))
# Base.show(io::IO, d::PowerMeasure{M,1}) where {M} = print(io, d.parent, " ^ ", length(d.xs))

# gentype(d::PowerMeasure{M,N}) where {M,N} = @inbounds Array{gentype(first(marginals(d))), N}

params(d::PowerMeasure) = params(first(marginals(d)))

# basemeasure(μ::PowerMeasure) = @inbounds basemeasure(first(μ.data))^size(μ.data)

@inline function basemeasure(d::PowerMeasure)
    basemeasure(d.parent)^d.axes
end

# Numeric flat variates of powers of measures with fixed stream sizes but
# no variate rank (e.g. tuple products) are streams:
@inline function mspace_ndims(::Type{<:PowerMeasure{M,A}}) where {M,A<:Tuple}
    _pwr_ndims(mspace_ndims(M), fieldcount(A), fixed_stream_size(M))
end
@inline _pwr_ndims(n::Integer, k::Integer, ::Any) = n + k
@inline _pwr_ndims(::NoMSpaceElementSize, ::Integer, ::True) = 1
@inline _pwr_ndims(n::NoMSpaceElementSize, ::Integer, ::False) = n

# Local measures of powers at nested variates are products of the local
# measures of the elements:
@inline localmeasure(μ::PowerMeasure, ::AbstractArray{<:Number}) = μ
function localmeasure(μ::PowerMeasure, x::AbstractArray)
    size(x) == _dynamic_dims(pwr_size(μ)) || return μ
    productmeasure(map(Base.Fix1(localmeasure, pwr_base(μ)), x))
end
@inline fixed_stream_size(::Type{<:PowerMeasure{M}}) where {M} = fixed_stream_size(M)

# The innermost base measure of nested powers and the total number of power
# dimensions:
@inline _pwr_unwrap(μ) = (μ, static(0))
@inline function _pwr_unwrap(μ::PowerMeasure)
    ν, n = _pwr_unwrap(pwr_base(μ))
    ν, n + static(length(pwr_axes(μ)))
end

# Batched kernels: the base kernel runs over the flat batch, the power then
# sums the leading dimensions of the result that belong to its axes.
@inline function _powered_kernel(f::F, μ::PowerMeasure, X) where {F}
    _check_pwr_batch(X, μ)
    _powered_kernel_impl(f, μ, X, _static_ndims(pwr_base(μ)))
end
@inline function _powered_kernel_impl(f::F, μ::PowerMeasure, X, ::Any) where {F}
    _sum_leading_dims(_batched_kernel(f, pwr_base(μ), X), static(length(pwr_axes(μ))))
end
# Numeric batches of powers of bases without a variate rank are batches of
# streams:
@inline function _powered_kernel_impl(::typeof(logdensityof_impl), μ::PowerMeasure, X::AbstractArray{<:Number}, ::NoMSpaceElementSize)
    _powered_stream_kernel(μ, X, fixed_stream_size(pwr_base(μ)))
end
function _powered_stream_kernel(μ::PowerMeasure, X::AbstractArray, ::True)
    ℓ, X_rest = batched_logdensityof_with_rest(μ, X, ())
    size(X_rest, 1) == 0 || _throw_stream_too_long()
    return ℓ
end
_powered_stream_kernel(μ::PowerMeasure, X::AbstractArray, ::False) = _streamwise_ld(logdensityof_impl, μ, X)

# Flat batches of powers have the power dimensions after the variate
# dimensions of the base measure (where the rank of the base is known):
@inline _check_pwr_batch(X::AbstractArray, μ::PowerMeasure) = _check_pwr_dims(X, _static_ndims(pwr_base(μ)), _dynamic_dims(pwr_size(μ)), false)
@inline _check_pwr_batch(::Any, ::PowerMeasure) = nothing
@inline function _check_pwr_dims(X::AbstractArray, ::StaticInteger{K}, dims::Dims, exact::Bool) where {K}
    n = length(dims)
    if (exact ? ndims(X) != K + n : ndims(X) < K + n) || ntuple(i -> size(X, K + i), Val(length(dims))) != dims
        _throw_size_mismatch()
    end
    return nothing
end
@inline _check_pwr_dims(::AbstractArray, ::NoMSpaceElementSize, ::Dims, ::Bool) = nothing
@inline _dynamic_dims(sz::SizeLike) = map(dynamic, _size_dims(sz))
@inline batched_logdensityof_impl(μ::PowerMeasure, X) = _powered_kernel(logdensityof_impl, μ, X)
@inline batched_logdensity_def(μ::PowerMeasure, X) = _powered_kernel(logdensity_def, μ, X)

# Point evaluation: flat variates are batches with zero batch dimensions,
# nested variates without flat storage sum the point densities of the base.
@inline _point_ld(f::F, μ::PowerMeasure, x::AbstractArray{<:Number}) where {F} = f(μ, x)
@inline logdensityof_impl(μ::PowerMeasure, x) = _powered_point(logdensityof_impl, μ, x)
@inline logdensity_def(μ::PowerMeasure, x) = _powered_point(logdensity_def, μ, x)

@inline function _powered_point(f::F, μ::PowerMeasure, x::AbstractArray{<:Number}) where {F}
    _point_result(_materialize(_batched_kernel(f, μ, x)), μ)
end
@inline function _powered_point(f::F, μ::PowerMeasure, x::AbstractArray) where {F}
    _check_pwr_shape(μ, x)
    _powered_point_nested(f, μ, x, _flat_storage(x))
end
@inline function _powered_point_nested(f::F, μ::PowerMeasure, x, x_flat::Union{AbstractArray,Tuple,NamedTuple}) where {F}
    _point_result(_materialize(_batched_kernel(f, μ, x_flat)), μ)
end
function _powered_point_nested(f::F, μ::PowerMeasure, x::AbstractArray, ::NoFlatStorage) where {F}
    ν = pwr_base(μ)
    sum(_PointLogd(f, ν), x; init = zero(_logd_numtype(x)))
end
@noinline function _powered_point(::F, ::PowerMeasure, x) where {F}
    throw(ArgumentError("Variates of powers of measures must be arrays"))
end

# Nested variates have the power's shape:
@inline function _check_pwr_shape(μ::PowerMeasure, x::AbstractArray)
    if maybestatic_size(x) != pwr_size(μ)
        _throw_size_mismatch()
    end
    return nothing
end

# Streams: a power consumes its size times the variates of the base measure
# and sums the base results over its axes. Bases without fixed variate
# sizes are consumed element by element, for single streams.
function batched_logdensityof_with_rest(μ::PowerMeasure, X::AbstractArray, sz::Dims)
    _powered_ld_with_rest(μ, X, sz, fixed_stream_size(pwr_base(μ)))
end
function batched_logdensityof_with_rest(μ::PowerMeasure, x::AbstractVector, sz::Tuple{})
    _powered_ld_with_rest(μ, x, sz, fixed_stream_size(pwr_base(μ)))
end
function _powered_ld_with_rest(μ::PowerMeasure, X::AbstractArray, sz::Dims, ::True)
    ℓ, X_rest = batched_logdensityof_with_rest(pwr_base(μ), X, (_dynamic_dims(pwr_size(μ))..., sz...))
    return _sum_leading_dims(ℓ, static(length(pwr_axes(μ)))), X_rest
end
function _powered_ld_with_rest(μ::PowerMeasure, x::AbstractVector, ::Tuple{}, ::False)
    ν = pwr_base(μ)
    ℓ = zero(_logd_numtype(x))
    x_rest = x
    for _ in 1:length(marginals(μ))
        ℓ_i, _, x_rest = logdensityof_with_rest(ν, x_rest)
        ℓ += ℓ_i
    end
    return ℓ, x_rest
end
@noinline function _powered_ld_with_rest(μ::PowerMeasure, ::AbstractArray, ::Dims, ::False)
    throw(ArgumentError("Batches of variate streams containing powers of measures of type $(nameof(typeof(pwr_base(μ)))) must be consumed stream by stream"))
end

# Support checks of powers run over the flat variate storage where the
# innermost base measure has scalar variates, elementwise otherwise:
@inline function insupport(μ::PowerMeasure, x::AbstractArray)
    ν, _ = _pwr_unwrap(μ)
    _powered_insupport(μ, x, _flat_storage(x), _static_ndims(ν))
end

@inline function _powered_insupport(μ::PowerMeasure, x, x_flat::AbstractArray, ::StaticInteger{0})
    ν, _ = _pwr_unwrap(μ)
    _all_insupport(broadcast(_insupport_bool ∘ Base.Fix1(insupport, ν), x_flat))
end
@inline function _powered_insupport(μ::PowerMeasure, x, x_flat::AbstractArray, ::StaticInteger{K}) where {K}
    ν, _ = _pwr_unwrap(μ)
    _all_insupport(map(_insupport_bool ∘ Base.Fix1(insupport, ν), sliced(x_flat, Val(K))))
end
@inline function _powered_insupport(μ::PowerMeasure, x, ::AbstractArray{<:Number}, ::NoMSpaceElementSize)
    NoFastInsupport{typeof(μ)}()
end
@inline _powered_insupport(μ::PowerMeasure, x, ::Any, ::Any) = _powered_insupport_elementwise(pwr_base(μ), x)

@inline function _powered_insupport_elementwise(ν, x::AbstractArray)
    _all_insupport(broadcast(_insupport_bool ∘ Base.Fix1(insupport, ν), x))
end

function insupport(μ::PowerMeasure, x)
    mapreduce(_insupport_bool ∘ Base.Fix1(insupport, pwr_base(μ)), _insupport_and, x)
end

@inline getdof(μ::PowerMeasure) = getdof(μ.parent) * size2length(axes2size(μ.axes))
@inline fast_dof(μ::PowerMeasure) = fast_dof(μ.parent) * size2length(axes2size(μ.axes))

# Static.SOneTo(0) is not static (yet):
@inline function getdof(::PowerMeasure{<:Any,<:NTuple{N,StaticOneToLike{0}}}) where {N}
    static(0)
end
@inline function fast_dof(::PowerMeasure{<:Any,<:NTuple{N,StaticOneToLike{0}}}) where {N}
    static(0)
end

# Variates may be nested arrays of the power's shape or their flat storage:
@propagate_inbounds function checked_arg(μ::PowerMeasure, x::AbstractArray{<:Any})
    @boundscheck _check_pwr_variate(μ, x)
    return x
end

@inline function _check_pwr_variate(μ::PowerMeasure, x::AbstractArray)
    if maybestatic_size(x) != pwr_size(μ)
        _check_pwr_flat(x, _static_ndims(pwr_base(μ)), _dynamic_dims(pwr_size(μ)))
    end
    return nothing
end
@inline _check_pwr_flat(x::AbstractArray, k::StaticInteger, dims::Dims) = _check_pwr_dims(x, k, dims, true)
@inline _check_pwr_flat(::AbstractArray, ::NoMSpaceElementSize, ::Dims) = _throw_size_mismatch()

checked_arg(μ::PowerMeasure, x::Any) = _throw_size_mismatch()

massof(m::PowerMeasure) = massof(m.parent)^dynamic(size2length(pwr_size(m)))


# Transport: the standard variate of a power is the flat vector of the
# standard variates of its innermost base measure, in the order of the flat
# variate storage. Batches transport over the flat storage `(base variate
# dims..., power dims..., batch dims...)`.

function batched_transport_to_std(::Type{S}, μ::PowerMeasure, X::AbstractArray) where {S<:StdMeasure}
    _check_pwr_batch(X, μ)
    _pwr_batched_to_std(S, μ, X, _static_ndims(pwr_base(μ)))
end
function batched_transport_to_std(::Type{S}, μ::PowerMeasure, X::Union{Tuple,NamedTuple}) where {S<:StdMeasure}
    _pwr_batched_to_std(S, μ, X, nothing)
end
@inline function _pwr_batched_to_std(::Type{S}, μ::PowerMeasure, X, ::Any) where {S}
    ν, n = _pwr_unwrap(μ)
    _merge_leading_dims(batched_transport_to_std(S, ν, X), static(1) + n)
end
# Numeric batches of powers of bases without a variate rank are batches of
# streams:
function _pwr_batched_to_std(::Type{S}, μ::PowerMeasure, X::AbstractArray{<:Number}, ::NoMSpaceElementSize) where {S}
    Z, X_rest = batched_transport_to_std_with_rest(S, μ, X, ())
    size(X_rest, 1) == 0 || _throw_stream_too_long()
    return Z
end

function batched_transport_from_std(::Type{S}, μ::PowerMeasure, Z::AbstractArray) where {S<:StdMeasure}
    ν, _ = _pwr_unwrap(μ)
    dims = _pwr_dims(μ)
    n_rows = _batch_dims(Z)[1]
    dof_ν = _base_dof(n_rows, prod(dims))
    dof_ν * prod(dims) == n_rows || _throw_std_length_mismatch()
    batched_transport_from_std(S, ν, _reshape_batch(Z, (dof_ν, dims..., Base.tail(_batch_dims(Z))...)))
end

# Empty powers leave the degrees of freedom of the base undetermined:
@inline _base_dof(n_rows::IntegerLike, n_pwr::IntegerLike) = n_rows ÷ max(n_pwr, one(n_pwr))

# All power dimensions of nested powers, innermost first:
@inline _pwr_dims(μ::PowerMeasure) = (_pwr_dims(pwr_base(μ))..., _size_dims(pwr_size(μ))...)
@inline _pwr_dims(ν) = ()

# Point transport: flat variates are batches with zero batch dimensions,
# nested variates without flat storage transport element by element.
function transport_to_std(::Type{S}, μ::PowerMeasure, x::AbstractArray) where {S<:StdMeasure}
    _pwr_to_std(S, μ, x, _flat_storage(x))
end
@inline function _pwr_to_std(::Type{S}, μ::PowerMeasure, x, x_flat::AbstractArray) where {S}
    _single_std(batched_transport_to_std(S, μ, x_flat))
end
function _pwr_to_std(::Type{S}, μ::PowerMeasure, x::AbstractArray, ::NoFlatStorage) where {S}
    _check_pwr_shape(μ, x)
    _flat_std_of(map(Base.Fix1(_ToStd{S}(), pwr_base(μ)), x))
end

function transport_from_std(::Type{S}, μ::PowerMeasure, z::AbstractVector) where {S<:StdMeasure}
    _pwr_variate(μ, batched_transport_from_std(S, μ, z))
end

# Streams: a power consumes the variates of its base measure with its size
# as additional multiplicity. Bases without fixed variate sizes are
# consumed element by element, for single streams.
function batched_transport_to_std_with_rest(::Type{S}, μ::PowerMeasure, X::AbstractArray, sz::Dims) where {S<:StdMeasure}
    _pwr_to_std_with_rest(S, μ, X, sz, fixed_stream_size(pwr_base(μ)))
end
@inline function _pwr_to_std_with_rest(::Type{S}, μ::PowerMeasure, X::AbstractArray, sz::Dims, ::True) where {S}
    batched_transport_to_std_with_rest(S, pwr_base(μ), X, (_dynamic_dims(pwr_size(μ))..., sz...))
end
function _pwr_to_std_with_rest(::Type{S}, μ::PowerMeasure, x::AbstractVector, ::Tuple{}, ::False) where {S}
    z, _, x_rest = transport_to_std_with_rest(S, μ, x)
    return z, x_rest
end
@noinline function _pwr_to_std_with_rest(::Type{S}, μ::PowerMeasure, ::AbstractArray, ::Dims, ::False) where {S}
    throw(ArgumentError("Batches of variate streams containing powers of measures of type $(nameof(typeof(pwr_base(μ)))) must be consumed stream by stream"))
end

function transport_to_std_with_rest(::Type{S}, μ::PowerMeasure, x::AbstractVector) where {S<:StdMeasure}
    _pwr_point_to_std_with_rest(S, μ, x, fixed_stream_size(pwr_base(μ)))
end
function _pwr_point_to_std_with_rest(::Type{S}, μ::PowerMeasure, x::AbstractVector, ::True) where {S}
    _pwr_point_to_std_with_rest(S, μ, x, _static_ndims(pwr_base(μ)))
end
function _pwr_point_to_std_with_rest(::Type{S}, μ::PowerMeasure, x::AbstractVector, ::StaticInteger) where {S}
    x_μ, x_rest = _consume_from_stream(x, _stream_consume_size(μ))
    return _as_stdstream(transport_to_std(S, μ, x_μ)), x_μ, x_rest
end
# Bases without a variate rank (tuple products) consume streams via the
# batched protocol:
function _pwr_point_to_std_with_rest(::Type{S}, μ::PowerMeasure, x::AbstractVector, ::NoMSpaceElementSize) where {S}
    z, x_rest = _pwr_to_std_with_rest(S, μ, x, (), static(true))
    x_μ, _ = _split_after(x, maybestatic_length(x) - maybestatic_length(x_rest))
    return z, x_μ, x_rest
end
function _pwr_point_to_std_with_rest(::Type{S}, μ::PowerMeasure, x::AbstractVector, ::False) where {S}
    ν = pwr_base(μ)
    zs = Vector{Any}(undef, length(marginals(μ)))
    x_rest = x
    for i in eachindex(zs)
        zs[i], _, x_rest = transport_to_std_with_rest(S, ν, x_rest)
    end
    x_μ, _ = _split_after(x, maybestatic_length(x) - maybestatic_length(x_rest))
    return reduce(vcat, [z for z in zs]), x_μ, x_rest
end

# Powers of measures without fast degrees of freedom transport their
# elements sequentially:
function transport_from_std_with_rest(::Type{S}, μ::PowerMeasure, z::AbstractVector) where {S<:StdMeasure}
    _pwr_from_std_with_rest(S, μ, z, fast_dof(μ))
end
@inline _pwr_from_std_with_rest(::Type{S}, μ, z, n::IntegerLike) where {S} = _from_std_with_rest_bydof(S, μ, z, n)
function _pwr_from_std_with_rest(::Type{S}, μ, z, ::AbstractNoDOF) where {S}
    _marginals_from_std_with_rest(S, marginals(μ), z)
end

# The stream length of a power with a base of fixed stream length:
@inline function _fixed_stream_length(μ::PowerMeasure)
    _fixed_stream_length(pwr_base(μ)) * prod(_dynamic_dims(pwr_size(μ)))
end

# The nested variate layout of a power over its flat storage, batches of
# tuple variates become struct arrays:
@inline _pwr_variate(μ::PowerMeasure, A::AbstractArray) = _pwr_variate_impl(μ, A)
@inline _pwr_variate(μ::PowerMeasure, A::Union{Tuple,NamedTuple}) = _pwr_variate_impl(μ, A)
@inline _pwr_variate_impl(μ::PowerMeasure, A) = _pwr_nest(pwr_base(μ), _pwr_variate(pwr_base(μ), A))
@inline _pwr_variate(ν, A::AbstractArray) = _nest_leaf(A, _static_ndims(ν))
@inline function _pwr_variate(ν::ProductMeasure{<:Tuple}, X::Tuple)
    StructArray(map((m, Xi) -> _nest_leaf(Xi, _static_ndims(m)), marginals(ν), X))
end
@inline function _pwr_variate(ν::ProductMeasure{<:NamedTuple{names}}, X::NamedTuple{names}) where {names}
    StructArray(NamedTuple{names}(map((m, Xi) -> _nest_leaf(Xi, _static_ndims(m)), values(marginals(ν)), values(X))))
end
@inline _nest_leaf(A::AbstractArray, ::StaticInteger{0}) = A
@inline _nest_leaf(A::AbstractArray, ::NoMSpaceElementSize) = A
@inline _nest_leaf(A::AbstractArray{<:Any,N}, ::StaticInteger{N}) where {N} = A
@inline _nest_leaf(A::AbstractArray, ::StaticInteger{K}) where {K} = sliced(A, Val(K))
@inline _pwr_nest(ν::PowerMeasure, B::AbstractArray) = sliced(B, Val(length(pwr_axes(ν))))
@inline _pwr_nest(ν, B::AbstractArray) = B

Adapt.adapt_structure(to, μ::PowerMeasure) = PowerMeasure(Adapt.adapt(to, pwr_base(μ)), pwr_axes(μ))
