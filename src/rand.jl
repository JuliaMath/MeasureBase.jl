# Random variate generation is parameterized by a `GenContext` carrying the
# random number generator, the numerical precision and the compute unit.
# Batches of variates are generated in flat form on the compute unit, single
# variates of powers are drawn as one batch and reshaped.

"""
    rand([rng::AbstractRNG], [T::Type{<:AbstractFloat}], μ::AbstractMeasure)
    rand(ctx::GenContext, μ::AbstractMeasure)

Generate a random variate of `μ`.

The generative context `ctx` (see `HeterogeneousComputing.GenContext`)
determines the random number generator, the numerical precision (`Float64`
by default) and the compute unit that array-valued variates are generated
on. The variates of powers of measures are generated in one batch.

Measure types should specialize [`MeasureBase.rand_impl`](@ref) and
[`MeasureBase.batched_rand_impl`](@ref) instead of `rand`.
"""
Base.rand(ctx::GenContext, μ::AbstractMeasure) = rand_impl(ctx, μ)

Base.rand(μ::AbstractMeasure) = rand(GenContext{Float64}(), μ)
Base.rand(rng::AbstractRNG, μ::AbstractMeasure) = rand(GenContext{Float64}(rng), μ)
Base.rand(::Type{T}, μ::AbstractMeasure) where {T<:AbstractFloat} = rand(GenContext{T}(), μ)
Base.rand(rng::AbstractRNG, ::Type{T}, μ::AbstractMeasure) where {T<:AbstractFloat} = rand(GenContext{T}(rng), μ)

@inline Random.rand!(d::AbstractMeasure, args...) = rand!(Random.default_rng(), d, args...)


"""
    MeasureBase.rand_impl(ctx::GenContext, μ)

Generate one random variate of `μ` in the generative context `ctx`.

The default implementation draws a variate of the preferred standard
measure of `μ` and transports it to `μ`. Measure types with a more direct
way of generating variates specialize `rand_impl`, and should specialize
[`MeasureBase.batched_rand_impl`](@ref) as well where batches can be
generated in a more direct way, too.
"""
function rand_impl end

function rand_impl(ctx::GenContext, μ)
    _rand_via_std(ctx, μ, preferred_stdmeasure(μ), fast_dof(μ), mspace_flatsize(μ))
end

@inline function _rand_via_std(ctx::GenContext, μ, ::Type{S}, ::IntegerLike, ::Tuple{}) where {S<:StdMeasure}
    transport_from_std(S, μ, rand_impl(ctx, S()))
end
@inline function _rand_via_std(ctx::GenContext, μ, ::Type{S}, n::IntegerLike, ::Any) where {S<:StdMeasure}
    transport_from_std(S, μ, _rand_std(ctx, S, (dynamic(n),)))
end
@inline function _rand_via_std(ctx::GenContext, μ, ::Type{AnyStdMeasure}, n::IntegerLike, sz)
    _rand_via_std(ctx, μ, StdUniform, n, sz)
end
function _rand_via_std(::GenContext, μ, ::Any, ::Any, ::Any)
    throw(ArgumentError("Random variate generation is not implemented for measures of type $(nameof(typeof(μ))), define MeasureBase.rand_impl"))
end


"""
    MeasureBase.batched_rand_impl(ctx::GenContext, μ, sz::Dims)

Generate a batch of random variates of `μ` of batch size `sz` in flat
form, an array of size `(flat variate dims..., sz...)` (see
[`MeasureBase.mspace_flatsize`](@ref)).

The default implementation draws a batch of variates of the preferred
standard measure of `μ` and transports it to `μ`, or generates the
variates one by one if `μ` has no standard transport.
"""
function batched_rand_impl end

function batched_rand_impl(ctx::GenContext, μ, sz::Dims)
    _batched_rand_via_std(ctx, μ, sz, preferred_stdmeasure(μ), fast_dof(μ), mspace_flatsize(μ))
end

function _batched_rand_via_std(ctx::GenContext, μ, sz::Dims, ::Type{S}, n::IntegerLike, ::SizeLike) where {S<:StdMeasure}
    batched_transport_from_std(S, μ, _rand_std(ctx, S, (dynamic(n), sz...)))
end
@inline function _batched_rand_via_std(ctx::GenContext, μ, sz::Dims, ::Type{AnyStdMeasure}, n::IntegerLike, sz_flat::SizeLike)
    _batched_rand_via_std(ctx, μ, sz, StdUniform, n, sz_flat)
end
@inline function _batched_rand_via_std(ctx::GenContext, μ, sz::Dims, ::Type{AnyStdMeasure}, n::IntegerLike, sz_flat::NoMSpaceElementSize)
    _batched_rand_via_std(ctx, μ, sz, StdUniform, n, sz_flat)
end
function _batched_rand_via_std(ctx::GenContext, μ, sz::Dims, ::Any, ::Any, ::SizeLike)
    _batched_rand_pointwise(ctx, μ, sz)
end
function _batched_rand_via_std(::GenContext, μ, ::Dims, ::Any, ::Any, ::NoMSpaceElementSize)
    throw(ArgumentError("Batched random variate generation requires measures of type $(nameof(typeof(μ))) to have a known variate size"))
end
function _batched_rand_via_std(::GenContext, μ, ::Dims, ::Type{S}, ::IntegerLike, ::NoMSpaceElementSize) where {S<:StdMeasure}
    throw(ArgumentError("Batched random variate generation requires measures of type $(nameof(typeof(μ))) to have a known variate size"))
end

function _batched_rand_pointwise(ctx::GenContext, μ, sz::Dims)
    _stack_variates(map(_ -> rand_impl(ctx, μ), CartesianIndices(sz)))
end

@inline _stack_variates(xs::AbstractArray{<:Number}) = xs
@inline _stack_variates(xs::AbstractArray{<:AbstractArray}) = stacked(xs)


# Bulk draws of standard variates on the compute unit:

@inline _rand_std(ctx::GenContext, ::Type{S}, dims::Dims) where {S<:StdMeasure} = batched_rand_impl(ctx, S(), dims)

@inline _rand_bulk(ctx::GenContext, sz::Dims) = rand(ctx, sz)
@inline _randn_bulk(ctx::GenContext, sz::Dims) = randn(ctx, sz)
@inline _randexp_bulk(ctx::GenContext, sz::Dims) = _randexp_bulk(ctx, sz, get_compute_unit(ctx))
@inline _randexp_bulk(ctx::GenContext, sz::Dims, ::CPUnit) = randexp(ctx, sz)
# Not all compute units provide exponential draws, derive them from uniform draws then:
@inline _randexp_bulk(ctx::GenContext, sz::Dims, ::AbstractComputeUnit) = -log1p.(-_rand_bulk(ctx, sz))

# Test values use a constant RNG, which only draws single values:
const _ConstantContext = GenContext{<:Any,<:Any,ConstantRNG}
@inline _rand_bulk(ctx::_ConstantContext, sz::Dims) = _const_bulk(ctx, rand(ConstantRNG(), get_precision(ctx)), sz)
@inline _randn_bulk(ctx::_ConstantContext, sz::Dims) = _const_bulk(ctx, randn(ConstantRNG(), get_precision(ctx)), sz)
@inline _randexp_bulk(ctx::_ConstantContext, sz::Dims) = _const_bulk(ctx, randexp(ConstantRNG(), get_precision(ctx)), sz)
@inline _const_bulk(ctx::GenContext, x, sz::Dims) = fill!(allocate_array(ctx, typeof(x), sz), x)

# A batch of copies of a constant variate:
function _const_batch(ctx::GenContext, x, sz::Dims)
    X = allocate_array(ctx, eltype(x), (size(x)..., sz...))
    X .= x
    return X
end
function _const_batch(ctx::GenContext, x::Number, sz::Dims)
    fill!(allocate_array(ctx, typeof(x), sz), x)
end
