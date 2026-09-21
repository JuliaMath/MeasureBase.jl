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

@inline Random.rand!(d::AbstractMeasure, args...) = rand!(Random.default_rng(), d, args...)


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
[`MeasureBase.batched_transport_from_std`](@ref)), or generates the
variates one by one via [`MeasureBase.rand_impl`](@ref) if `μ` has no
standard transport.
"""
function batched_rand_impl end

"""
    MeasureBase.rand_impl(ctx::GenContext, μ)

Generate one random variate of `μ` in the generative context `ctx`.

The default implementation generates a batch with zero batch dimensions
via [`MeasureBase.batched_rand_impl`](@ref). Measure types with a more
direct way of generating single variates may specialize `rand_impl`.
"""
function rand_impl end

# The marker tells the defaults whether `rand_impl` may be specialized
# for the measure (coming from the default `rand_impl` itself, it is not):
struct _NoRandImpl end
struct _MaybeRandImpl end

@inline rand_impl(ctx::GenContext, μ) = _rand_default(ctx, μ, (), _NoRandImpl())
@inline batched_rand_impl(ctx::GenContext, μ, sz::SizeLike) = _rand_default(ctx, μ, sz, _MaybeRandImpl())

@inline _rand_default(ctx::GenContext, μ, sz::SizeLike, m) = _rand_via_std(ctx, μ, sz, preferred_stdmeasure(μ), m)

@inline function _rand_via_std(ctx::GenContext, μ, sz::SizeLike, ::Type{S}, m) where {S<:StdMeasure}
    _rand_via_std_dof(ctx, μ, sz, S, fast_dof(μ), m)
end
@inline _rand_via_std(ctx::GenContext, μ, sz::SizeLike, ::Type{AnyStdMeasure}, m) = _rand_via_std(ctx, μ, sz, StdUniform, m)
@inline _rand_via_std(ctx::GenContext, μ, sz::SizeLike, ::Any, m) = _rand_pointwise(ctx, μ, sz, m)

function _rand_via_std_dof(ctx::GenContext, μ, sz::SizeLike, ::Type{S}, n::IntegerLike, ::Any) where {S<:StdMeasure}
    convert_realtype(get_precision(ctx), batched_transport_from_std(S, μ, _rand_std(ctx, S, (n, size_dims(sz)...))))
end
@inline _rand_via_std_dof(ctx::GenContext, μ, sz::SizeLike, ::Type, ::Any, m) = _rand_pointwise(ctx, μ, sz, m)

# Variates generated one by one, stacked into a flat batch:
@inline _rand_pointwise(ctx::GenContext, μ, sz::SizeLike, ::Any) = _batched_rand_pointwise(ctx, μ, sz)
@inline _rand_pointwise(ctx::GenContext, μ, ::Tuple{}, ::_MaybeRandImpl) = rand_impl(ctx, μ)
@noinline function _rand_pointwise(::GenContext, μ, ::Tuple{}, ::_NoRandImpl)
    throw(ArgumentError("Random variate generation is not implemented for measures of type $(nameof(typeof(μ))), define MeasureBase.batched_rand_impl or MeasureBase.rand_impl"))
end

function _batched_rand_pointwise(ctx::GenContext, μ, sz::SizeLike)
    _stack_variates(map(_ -> rand_impl(ctx, μ), CartesianIndices(asnonstatic(sz))))
end
@inline _batched_rand_pointwise(ctx::GenContext, μ, ::Tuple{}) = rand_impl(ctx, μ)

@inline _stack_variates(xs::AbstractArray{<:Number}) = xs
@inline _stack_variates(xs::AbstractArray{<:AbstractArray}) = stacked(map(_stack_variates, xs))
@inline _stack_variates(xs::AbstractArray{<:Union{Tuple,NamedTuple}}) = StructArrays.components(StructArray(xs))


# Bulk draws of standard variates on the compute unit, single draws for
# zero batch dimensions:

@inline _rand_std(ctx::GenContext, ::Type{S}, dims::SizeLike) where {S<:StdMeasure} = batched_rand_impl(ctx, S(), dims)

@inline _rand_bulk(ctx::GenContext, sz::SizeLike) = _bulk_draw(rand, ctx, sz)
@inline _randn_bulk(ctx::GenContext, sz::SizeLike) = _bulk_draw(randn, ctx, sz)
@inline _randexp_bulk(ctx::GenContext, sz::SizeLike) = _randexp_bulk(ctx, sz, get_compute_unit(ctx))
@inline _randexp_bulk(ctx::GenContext, sz::SizeLike, ::CPUnit) = _bulk_draw(randexp, ctx, sz)
# Not all compute units provide exponential draws, derive them from uniform draws then:
@inline _randexp_bulk(ctx::GenContext, sz::SizeLike, ::AbstractComputeUnit) = -log1p.(-_rand_bulk(ctx, sz))

# Fully static batch sizes draw static arrays on the CPU, so that variates
# of statically sized measures are allocation-free. Other compute units
# allocate their own arrays.
@inline _bulk_draw(f::F, ctx::GenContext, sz::SizeLike) where {F} = f(ctx, asnonstatic(sz))
@inline _bulk_draw(f::F, ctx::GenContext, sz::StaticSizeLike) where {F} =
    _bulk_draw(f, ctx, sz, get_compute_unit(ctx))
@inline _bulk_draw(f::F, ctx::GenContext, sz::StaticSizeLike, ::AbstractComputeUnit) where {F} =
    f(ctx, asnonstatic(sz))
@inline _bulk_draw(f::F, ctx::GenContext, sz::StaticSizeLike, ::CPUnit) where {F} =
    f(get_rng(ctx), staticarray_type(get_precision(ctx), canonical_size(sz)))

@inline _rand_bulk(ctx::GenContext, ::Tuple{}) = rand(get_rng(ctx), get_precision(ctx))
@inline _randn_bulk(ctx::GenContext, ::Tuple{}) = randn(get_rng(ctx), get_precision(ctx))
@inline _randexp_bulk(ctx::GenContext, ::Tuple{}) = randexp(get_rng(ctx), get_precision(ctx))

# Test values use a constant RNG, which only draws single values:
const _ConstantContext = GenContext{<:AbstractFloat,<:AbstractComputeUnit,ConstantRNG}
@inline _rand_bulk(ctx::_ConstantContext, sz::SizeLike) = _const_bulk(ctx, rand(ConstantRNG(), get_precision(ctx)), sz)
@inline _randn_bulk(ctx::_ConstantContext, sz::SizeLike) = _const_bulk(ctx, randn(ConstantRNG(), get_precision(ctx)), sz)
@inline _randexp_bulk(ctx::_ConstantContext, sz::SizeLike) = _const_bulk(ctx, randexp(ConstantRNG(), get_precision(ctx)), sz)
@inline _rand_bulk(ctx::_ConstantContext, ::Tuple{}) = rand(ConstantRNG(), get_precision(ctx))
@inline _randn_bulk(ctx::_ConstantContext, ::Tuple{}) = randn(ConstantRNG(), get_precision(ctx))
@inline _randexp_bulk(ctx::_ConstantContext, ::Tuple{}) = randexp(ConstantRNG(), get_precision(ctx))
@inline _const_bulk(ctx::GenContext, x, sz::SizeLike) = _const_bulk(ctx, x, sz, get_compute_unit(ctx))
@inline _const_bulk(ctx::GenContext, x, sz::SizeLike, ::AbstractComputeUnit) =
    fill!(allocate_array(ctx, typeof(x), asnonstatic(sz)), x)
@inline _const_bulk(ctx::GenContext, x, sz::StaticSizeLike, ::CPUnit) = maybestatic_fill(x, sz)

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
