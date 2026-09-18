# Batched transport over flat batches of variates: the leading dimensions of
# a batch are the variate dimensions (see `mspace_flatsize`), all further
# dimensions are batch dimensions. Streams of standard variates are batches
# `(dof, batch dims...)`, consumed along their first dimension.

"""
    MeasureBase.batched_transport_to_std(::Type{S}, μ, X::AbstractArray)

Batched form of [`MeasureBase.transport_to_std`](@ref): transports the
flat batch `X` of variates of `μ` to a batch `(getdof(μ), batch dims...)`
of variates of the standard measure type `S`.

The default implementation broadcasts the point transport for measures
with scalar variates and maps it over the variate slices of `X`
otherwise, in a host loop.
"""
function batched_transport_to_std end

function batched_transport_to_std(::Type{S}, μ, X::AbstractArray) where {S<:StdMeasure}
    _batched_to_std(S, μ, X, mspace_flatsize(μ))
end

@inline function _batched_to_std(::Type{S}, μ, X::AbstractArray, ::Tuple{}) where {S}
    _as_stdstream_batch(broadcast(Base.Fix1(_ToStd{S}(), μ), X))
end

function _batched_to_std(::Type{S}, μ, X::AbstractArray, sz::SizeLike) where {S}
    _check_flatsize(X, sz)
    stacked(map(Base.Fix1(_ToStd{S}(), μ), sliced(X, Val(length(sz)))))
end

function _batched_to_std(::Type{S}, μ, ::AbstractArray, ::NoMSpaceElementSize) where {S}
    throw(ArgumentError("Batched transport requires measures of type $(nameof(typeof(μ))) to have a known variate size"))
end

# Standard variates of scalar-variate measures form the first dimension:
@inline _as_stdstream_batch(Z::AbstractArray) = reshape(Z, (1, size(Z)...))
@inline _drop_stdstream_dim(Z::AbstractArray) = reshape(Z, Base.tail(size(Z)))


"""
    MeasureBase.batched_transport_from_std(::Type{S}, μ, Z::AbstractArray)

Batched form of [`MeasureBase.transport_from_std`](@ref): transports the
batch `Z` of variates of the standard measure type `S`, of size
`(getdof(μ), batch dims...)`, to a flat batch of variates of `μ`.

The default implementation broadcasts the point transport for measures
with scalar variates and maps it over the columns of `Z` otherwise, in a
host loop.
"""
function batched_transport_from_std end

function batched_transport_from_std(::Type{S}, μ, Z::AbstractArray) where {S<:StdMeasure}
    _batched_from_std(S, μ, Z, mspace_flatsize(μ))
end

@inline function _batched_from_std(::Type{S}, μ, Z::AbstractArray, ::Tuple{}) where {S}
    broadcast(Base.Fix1(_FromStd{S}(), μ), _drop_stdstream_dim(Z))
end

function _batched_from_std(::Type{S}, μ, Z::AbstractArray, sz::SizeLike) where {S}
    xs = map(Base.Fix1(_FromStd{S}(), μ), sliced(Z, Val(1)))
    reshape(stacked(xs), (map(dynamic, _size_dims(sz))..., Base.tail(size(Z))...))
end

function _batched_from_std(::Type{S}, μ, ::AbstractArray, ::NoMSpaceElementSize) where {S}
    throw(ArgumentError("Batched transport requires measures of type $(nameof(typeof(μ))) to have a known variate size"))
end


"""
    MeasureBase.batched_transport_to_std_with_rest(::Type{S}, μ, X::AbstractArray)

Batched form of [`MeasureBase.transport_to_std_with_rest`](@ref) for a
batch `X` of flat vector streams (first dimension along the streams).

Returns a tuple `(Z, X_μ, X_rest)` of the batch of standard variates, the
batch of variates of `μ` consumed from the streams and the unconsumed rest
of the streams.
"""
function batched_transport_to_std_with_rest end

function batched_transport_to_std_with_rest(::Type{S}, μ, X::AbstractArray) where {S<:StdMeasure}
    vsz = _stream_consume_size(μ)
    X_μ, X_rest = _batched_consume(X, vsz, ())
    X_v = _consumed_variates(X_μ, vsz)
    return batched_transport_to_std(S, μ, X_v), X_v, X_rest
end
@inline _consumed_variates(X_μ::AbstractArray, ::Tuple{}) = _drop_stdstream_dim(X_μ)
@inline _consumed_variates(X_μ::AbstractArray, ::SizeLike) = X_μ


"""
    MeasureBase.batched_transport_from_std_with_rest(::Type{S}, μ, Z::AbstractArray)

Batched form of [`MeasureBase.transport_from_std_with_rest`](@ref) for a
batch `Z` of streams of standard variates (first dimension along the
streams).

Returns a tuple `(X, Z_rest)` of the flat batch of variates of `μ` and the
unconsumed rest of the streams.
"""
function batched_transport_from_std_with_rest end

function batched_transport_from_std_with_rest(::Type{S}, μ, Z::AbstractArray) where {S<:StdMeasure}
    _batched_from_std_with_rest_bydof(S, μ, Z, fast_dof(μ))
end

function _batched_from_std_with_rest_bydof(::Type{S}, μ, Z::AbstractArray, n::IntegerLike) where {S}
    Z_μ, Z_rest = _batched_split(Z, dynamic(n))
    return batched_transport_from_std(S, μ, Z_μ), Z_rest
end

function _batched_from_std_with_rest_bydof(::Type{S}, μ, ::AbstractArray, ::AbstractNoDOF) where {S}
    throw(ArgumentError("Batched transport from standard measures requires measures of type $(nameof(typeof(μ))) to implement MeasureBase.batched_transport_from_std_with_rest"))
end

# A flat batch of variates as a batch of streams, the variate dimensions
# merged into the first dimension:
@inline function _as_stream_batch(X::AbstractArray, sz::SizeLike)
    n = length(sz)
    batch_dims = ntuple(i -> size(X, n + i), Val(ndims(X) - n))
    reshape(X, (dynamic(size2length(sz)), batch_dims...))
end


"""
    MeasureBase.batched_transport_def(ν, μ, X::AbstractArray)

Transport the flat batch `X` of variates of `μ` to a flat batch of
variates of `ν`, via the standard measure type the preferences of `ν` and
`μ` promote to. Specialize for pairs of measure types with a direct
batched transport.
"""
function batched_transport_def end

function batched_transport_def(ν, μ, X::AbstractArray)
    S = _transport_pivot(ν, μ)
    Z = batched_transport_to_std(S, μ, X)
    Y, Z_rest = batched_transport_from_std_with_rest(S, ν, Z)
    if size(Z_rest, 1) != 0
        throw(ArgumentError("Degrees of freedom of source and target measure of a transport don't match"))
    end
    return Y
end


# Broadcasting a transport function over an array of variates with flat
# storage, or over the flat storage of a batch, transports the batch as a
# whole. Fused broadcast arguments are materialized first, static arrays
# are transported point by point.
function Broadcast.broadcasted(f::TransportFunction, X::AbstractArray)
    _broadcast_transport(f, X, _flat_storage(X), mspace_flatsize(f.μ), mspace_flatsize(f.ν))
end

function Broadcast.broadcasted(f::TransportFunction, bc::Broadcast.Broadcasted)
    Broadcast.broadcasted(f, Broadcast.materialize(bc))
end

Broadcast.broadcasted(f::TransportFunction, X::StaticArray) = map(_Pointwise(f), X)

function _broadcast_transport(f::TransportFunction, X, X_flat::AbstractArray, sz_μ::SizeLike, sz_ν::SizeLike)
    _check_flatsize(X_flat, sz_μ)
    Y_flat = batched_transport_def(f.ν, f.μ, X_flat)
    return _batch_variates(Y_flat, f.ν)
end

_broadcast_transport(f::TransportFunction, X, ::Any, ::Any, ::Any) = map(_Pointwise(f), X)

# Prevents re-entering the broadcast hook from `map` implementations that
# broadcast (e.g. GPU arrays):
struct _Pointwise{F} <: Function
    f::F
end
@inline (p::_Pointwise)(x) = p.f(x)

# The batch of variates in the layout of the target measure over the flat
# result, nested powers included:
@inline _batch_variates(Y::AbstractArray, ν) = _nest_leaf(Y, mspace_flatsize(ν))
@inline _batch_variates(Y::AbstractArray, ν::PowerMeasure) = sliced(_pwr_variate(ν, Y), Val(length(pwr_axes(ν))))
