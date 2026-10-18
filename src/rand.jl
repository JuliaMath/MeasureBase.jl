# Random variate generation is parameterized by a `GenContext` carrying the
# random number generator, the numerical precision and the compute unit.
# Variates are generated as flat batches `(variate dims..., batch dims...)`
# on the compute unit, a single variate is a batch with zero batch
# dimensions.

"""
    rand([rng::AbstractRNG], [T::Type{<:AbstractFloat}], μ::AbstractMeasure)
    rand(ctx::GenContext, μ::AbstractMeasure)

Generate a random variate of `μ`.

The generative context `ctx` (see `HeterogeneousComputing.GenContext`)
determines the random number generator, the numerical precision (`Float64`
by default) and the compute unit that array-valued variates are generated
on. Variates of powers of measures are generated as one flat batch of
variates of the base measure.

Measure types should specialize [`MeasureBase.batched_rand_impl`](@ref)
instead of `rand`.
"""
Base.rand(ctx::GenContext, μ::AbstractMeasure) = rand_impl(ctx, μ)

Base.rand(μ::AbstractMeasure) = rand(GenContext{Float64}(), μ)
Base.rand(rng::AbstractRNG, μ::AbstractMeasure) = rand(GenContext{Float64}(rng), μ)
Base.rand(::Type{T}, μ::AbstractMeasure) where {T<:AbstractFloat} = rand(GenContext{T}(), μ)
Base.rand(rng::AbstractRNG, ::Type{T}, μ::AbstractMeasure) where {T<:AbstractFloat} = rand(GenContext{T}(rng), μ)


"""
    MeasureBase.batched_rand_impl(ctx::GenContext, μ, sz::SizeLike)

Generate a batch of random variates of `μ` of batch size `sz` in flat
form, an array `(variate dims..., sz...)`, or a single variate for
`sz == ()`. Batches of tuple and named tuple variates are tuples resp.
named tuples of batches. Fully static batch sizes give static arrays on
the CPU.

This is the primary extension point for random variate generation. The
default implementation draws a batch of variates of the preferred
standard measure of `μ` and transports it to `μ` (see
[`MeasureBase.batched_transport_from_std`](@ref)). Measures without a
standard transport get batches generated variate by variate via
[`MeasureBase.rand_impl`](@ref).
"""
function batched_rand_impl end

"""
    MeasureBase.rand_impl(ctx::GenContext, μ)

Generate one random variate of `μ` in the generative context `ctx`.

The default implementation calls
[`batched_rand_impl(ctx, μ, ())`](@ref MeasureBase.batched_rand_impl).
Measure types may specialize `rand_impl` if they have a more direct way of
generating single variates, measure types without a standard transport
may specialize `rand_impl` instead of `batched_rand_impl`. Structural
measures generate single variates of their components via `rand_impl`.
"""
function rand_impl end

# The defaults call each other only for batches: single variates of
# measures without a standard transport must come from a specialization of
# either function.
@inline rand_impl(ctx::GenContext, μ) = batched_rand_impl(ctx, μ, ())
@inline batched_rand_impl(ctx::GenContext, μ, sz::SizeLike) = _rand_default(ctx, μ, sz)

@inline _rand_default(ctx::GenContext, μ, sz::SizeLike) = _rand_via_std(ctx, μ, sz, preferred_stdmeasure(μ))
# Transitional, for callers that still pass the former `rand_impl` marker:
struct _NoRandImpl end
@inline _rand_default(ctx::GenContext, μ, sz::SizeLike, ::_NoRandImpl) = _rand_default(ctx, μ, sz)

@inline function _rand_via_std(ctx::GenContext, μ, sz::SizeLike, ::Type{S}) where {S<:StdMeasure}
    _rand_via_std_dof(ctx, μ, sz, S, fast_dof(μ))
end
@inline _rand_via_std(ctx::GenContext, μ, sz::SizeLike, ::Type{AnyStdMeasure}) = _rand_via_std(ctx, μ, sz, StdUniform)
@inline _rand_via_std(ctx::GenContext, μ, sz::SizeLike, ::Any) = _rand_pointwise(ctx, μ, sz)

function _rand_via_std_dof(ctx::GenContext, μ, sz::SizeLike, ::Type{S}, n::IntegerLike) where {S<:StdMeasure}
    convert_realtype(get_precision(ctx), _transported_std_draw(ctx, μ, sz, S, n, _static_ndims(μ)))
end
@inline _rand_via_std_dof(ctx::GenContext, μ, sz::SizeLike, ::Type, ::Any) = _rand_pointwise(ctx, μ, sz)

@inline _transported_std_draw(ctx::GenContext, μ, sz::SizeLike, ::Type{S}, n::IntegerLike, ::Any) where {S} =
    batched_transport_from_std(S, μ, _rand_std(ctx, S, (n, size_tuple(sz)...)))
# Single scalar variates from single standard variates, device arrays of
# standard variates can't be indexed on the host:
@inline _transported_std_draw(ctx::GenContext, μ, ::Tuple{}, ::Type{S}, ::IntegerLike, ::StaticInteger{0}) where {S} =
    transport_from_std(S, μ, _rand_std(ctx, S, ()))

# Variates generated one by one, stacked into a flat batch:
@inline _rand_pointwise(ctx::GenContext, μ, sz::SizeLike) = _batched_rand_pointwise(ctx, μ, sz)
@noinline _rand_pointwise(::GenContext, μ, ::Tuple{}) = _throw_no_rand_impl(μ)
@noinline function _throw_no_rand_impl(μ)
    throw(ArgumentError("Random variate generation is not implemented for measures of type $(nameof(typeof(μ))), define MeasureBase.batched_rand_impl or MeasureBase.rand_impl"))
end

function _batched_rand_pointwise(ctx::GenContext, μ, sz::SizeLike)
    _stack_variates(map(_ -> rand_impl(ctx, μ), CartesianIndices(asnonstatic(sz))), μ)
end
@inline _batched_rand_pointwise(ctx::GenContext, μ, ::Tuple{}) = rand_impl(ctx, μ)

@inline _stack_variates(xs::AbstractArray, μ) = _stack_variates(xs)
# Empty draws of measures without a sampler never call it, so nothing tells the variate type:
@noinline _stack_variates(::AbstractArray{Union{}}, μ) = _throw_no_rand_impl(μ)
@inline _stack_variates(xs::AbstractArray{<:Number}) = xs
@inline _stack_variates(xs::AbstractArray{<:AbstractArray}) = stacked(map(_stack_variates, xs))
@inline _stack_variates(xs::AbstractArray{<:Union{Tuple,NamedTuple}}) = StructArrays.components(StructArray(xs))


# Bulk draws of standard variates on the compute unit, single draws for
# zero batch dimensions:

@inline _rand_std(ctx::GenContext, ::Type{S}, dims::SizeLike) where {S<:StdMeasure} = batched_rand_impl(ctx, S(), dims)

@inline _rand_bulk(ctx::GenContext, sz::SizeLike) = _bulk_draw(rand, ctx, sz)
@inline _randn_bulk(ctx::GenContext, sz::SizeLike) = _bulk_draw(randn, ctx, sz)
@inline _randexp_bulk(ctx::GenContext, sz::SizeLike) = _bulk_draw(randexp, ctx, sz)

@inline _bulk_draw(f::F, ctx::GenContext, ::Tuple{}) where {F} = f(ctx)
@inline _bulk_draw(f::F, ctx::GenContext, sz::SizeLike) where {F} = _bulk_draw(f, ctx, sz, get_rng(ctx))
@inline _bulk_draw(f::F, ctx::GenContext, sz::SizeLike, ::AbstractRNG) where {F} = _maybestatic_draw(f, ctx, sz)
# Test values use a constant RNG, which only draws single values:
@inline _bulk_draw(f::F, ctx::GenContext, sz::SizeLike, ::ConstantRNG) where {F} = _const_bulk(ctx, f(ctx), sz)

# Fully static batch sizes draw static arrays on the host, so that variates
# of statically sized measures are allocation-free:
@inline _maybestatic_draw(f::F, ctx::GenContext, sz::SizeLike) where {F} = f(ctx, asnonstatic(sz))
@inline function _maybestatic_draw(f::F, ctx::GenContext, sz::StaticSizeLike) where {F}
    _on_host(ctx) ? f(get_rng(ctx), staticarray_type(get_precision(ctx), canonical_size(sz))) : f(ctx, asnonstatic(sz))
end

@inline _on_host(ctx::GenContext) = is_host_unit(get_compute_unit(ctx))

@inline _const_bulk(ctx::GenContext, x, sz::SizeLike) = fill_array(ctx, x, asnonstatic(sz))
@inline _const_bulk(ctx::GenContext, x, sz::StaticSizeLike) =
    _on_host(ctx) ? maybestatic_fill(x, sz) : fill_array(ctx, x, asnonstatic(sz))

# A mask over the batch dimensions, aligned with a flat batch of variates
# of rank `k`:
@inline _batch_mask(mask::Number, ::Any) = mask
@inline _batch_mask(mask::AbstractArray, ::StaticInteger{0}) = mask
@inline function _batch_mask(mask::AbstractArray, ::StaticInteger{K}) where {K}
    reshape(mask, (ntuple(_ -> 1, Val(K))..., size(mask)...))
end

# A batch of copies of a constant variate:
function _const_batch(ctx::GenContext, x, sz::SizeLike)
    X = allocate_array(ctx, eltype(x), (size(x)..., asnonstatic(sz)...))
    X .= x
    return X
end
@inline _const_batch(ctx::GenContext, x::Number, sz::SizeLike) = _const_bulk(ctx, x, sz)
@inline _const_batch(::GenContext, x::Number, ::Tuple{}) = x
