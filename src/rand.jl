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

Batches of variates are variates of powers of `μ`, `rand(μ^n)` resp.
`rand(μ^dims)`. Variates follow the normalized law of `μ`, weights don't
affect them.
Binds are sampled ancestrally from the normalized laws of their components,
see [`mbind`](@ref).

Floating-point variates have the precision of the context, also through
transformations. Point masses generate their atom exactly, and integer
variates (e.g. of counting measures) stay integers.

Measure types should specialize [`MeasureBase.batched_rand_impl`](@ref)
instead of `rand`.
"""
Base.rand(ctx::GenContext, μ::AbstractMeasure) = rand_impl(ctx, μ)

Base.rand(μ::AbstractMeasure) = rand(GenContext{Float64}(), μ)
Base.rand(rng::AbstractRNG, μ::AbstractMeasure) = rand(GenContext{Float64}(rng), μ)
Base.rand(::Type{T}, μ::AbstractMeasure) where {T<:AbstractFloat} = rand(GenContext{T}(), μ)
Base.rand(rng::AbstractRNG, ::Type{T}, μ::AbstractMeasure) where {T<:AbstractFloat} = rand(GenContext{T}(rng), μ)

# Batches of variates are variates of powers:
Base.rand(μ::AbstractMeasure, dims::Dims) = _throw_rand_dims()
Base.rand(μ::AbstractMeasure, dim::Integer, dims::Integer...) = _throw_rand_dims()
Base.rand(rng::AbstractRNG, μ::AbstractMeasure, dims::Dims) = _throw_rand_dims()
Base.rand(rng::AbstractRNG, μ::AbstractMeasure, dim::Integer, dims::Integer...) = _throw_rand_dims()
@noinline _throw_rand_dims() =
    throw(ArgumentError("rand(μ, dims...) is not supported for measures, use rand(μ^dims) for a batch of variates of μ"))


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
    _rand_via_std_dof(ctx, μ, sz, S, _split_dof(μ))
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
# A single static degree of freedom without a declared rank may be a
# scalar or a one-element array variate:
@noinline _transported_std_draw(::GenContext, μ, ::SizeLike, ::Type, ::StaticInteger{1}, ::NoMSpaceElementSize) =
    _throw_undeclared_rank("Random variate generation", μ)

# Variates generated one by one in a host loop, stacked into a flat batch
# and moved to the compute unit:
@inline _rand_pointwise(ctx::GenContext, μ, sz::SizeLike) = _batched_rand_pointwise(ctx, μ, sz)
@noinline _rand_pointwise(::GenContext, μ, ::Tuple{}) = _throw_no_rand_impl(μ)
@noinline function _throw_no_rand_impl(μ)
    throw(ArgumentError("Random variate generation is not implemented for measures of type $(nameof(typeof(μ))), define MeasureBase.batched_rand_impl or MeasureBase.rand_impl"))
end

function _batched_rand_pointwise(ctx::GenContext, μ, sz::SizeLike)
    xs = map(_ -> rand_impl(ctx, μ), CartesianIndices(asnonstatic(sz)))
    _on_compute_unit(ctx, _stack_variates(xs, μ))
end
@inline _batched_rand_pointwise(ctx::GenContext, μ, ::Tuple{}) = rand_impl(ctx, μ)

@inline _stack_variates(xs::AbstractArray, μ) = _stack_variates(xs)
# Empty draws of measures without a sampler never call it, so nothing tells the variate type:
@noinline _stack_variates(::AbstractArray{Union{}}, μ) = _throw_no_rand_impl(μ)
@inline _stack_variates(xs::AbstractArray{<:Number}) = xs
@inline _stack_variates(xs::AbstractArray{<:AbstractArray}) = stacked(map(_stack_variates, xs))
@inline _stack_variates(xs::AbstractArray{<:Union{Tuple,NamedTuple}}) = StructArrays.components(StructArray(xs))
# Empty batches of array variates keep the variate dimensions:
function _stack_variates(xs::AbstractArray{<:AbstractArray{T}}, μ) where {T<:Number}
    isempty(xs) ? similar(xs, T, (_empty_batch_flatsize(μ)..., size(xs)...)) : _stack_variates(xs)
end

# Numeric batches generated on the host, moved to the compute unit of the
# context. Fills and device arrays stay as they are, and so do traced
# numbers:
@inline _on_compute_unit(ctx::GenContext, X) = _on_host(ctx) ? X : _to_compute_unit(get_compute_unit(ctx), X)
@inline _to_compute_unit(unit, X::Union{Tuple,NamedTuple}) = map(Base.Fix1(_to_compute_unit, unit), X)
@inline _to_compute_unit(unit, X) = X
function _to_compute_unit(unit, X::Array{T}) where {T<:Number}
    isbitstype(T) || return X
    Y = allocate_array(unit, T, size(X))
    copyto!(Y, X)
    return Y
end


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

# Variates that only host code can generate are drawn in a host context
# derived from the context and moved to its compute unit. The RNG of the
# derived context is seeded with numbers drawn from the RNG of the context
# (see `adapt` for `GenContext`), host contexts stay as they are. Compiled
# code can't derive contexts:
@inline _host_context(ctx::GenContext) = Adapt.adapt(CPUnit(), ctx)
@inline _can_derive_host_context(::GenContext) = true
@inline _draws_in_host_context(ctx::GenContext) = !_on_host(ctx) && _can_derive_host_context(ctx)

# Marginal types whose variates, numbers or arrays of numbers, only host
# code can generate, one variate at a time:
@inline _host_drawn_marginals(::Type) = false

# Constant fills keep the type of the constant (point masses are exact):
@inline _const_bulk(ctx::GenContext, x, sz::SizeLike) = fill_array(get_compute_unit(ctx), x, asnonstatic(sz))
@inline _const_bulk(ctx::GenContext, x, sz::StaticSizeLike) =
    _on_host(ctx) ? maybestatic_fill(x, sz) : fill_array(get_compute_unit(ctx), x, asnonstatic(sz))

# A mask over the batch dimensions, aligned with a flat batch of variates
# of rank `k`:
@inline _batch_mask(mask::Number, ::Any) = mask
@inline _batch_mask(mask::AbstractArray, ::StaticInteger{0}) = mask
@inline function _batch_mask(mask::AbstractArray, ::StaticInteger{K}) where {K}
    reshape(mask, (ntuple(_ -> 1, Val(K))..., size(mask)...))
end

# A batch of copies of a constant variate. On the host, batches of numbers
# and fills are fills (static arrays for static sizes), so that constant
# components of batches stay recognizable:
@inline function _const_batch(ctx::GenContext, x, sz::SizeLike)
    _on_host(ctx) && _has_fill_batch(x) ? _fill_batch(x, sz) : _dense_const_batch(ctx, x, sz)
end
@inline _const_batch(::GenContext, x::Number, ::Tuple{}) = x

@inline _has_fill_batch(::Union{Number,FillArrays.AbstractFill}) = true
@inline _has_fill_batch(::Any) = false

@inline _fill_batch(x::Number, sz::SizeLike) = maybestatic_fill(x, sz)
@inline function _fill_batch(x::FillArrays.AbstractFill, sz::SizeLike)
    maybestatic_fill(FillArrays.getindex_value(x), (size(x)..., size_tuple(sz)...))
end

@inline _dense_const_batch(ctx::GenContext, x::Number, sz::SizeLike) = _const_bulk(ctx, x, sz)
function _dense_const_batch(ctx::GenContext, x, sz::SizeLike)
    X = allocate_array(ctx, eltype(x), (size(x)..., asnonstatic(sz)...))
    X .= x
    return X
end
