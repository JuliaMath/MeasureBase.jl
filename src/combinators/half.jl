export Half

"""
    Half(μ)

The distribution of `abs(x)` for `x` drawn from `μ`, for measures `μ` with
scalar variates that are symmetric about zero.
"""
struct Half{M} <: AbstractMeasure
    parent::M

    function Half(parent::M) where {M}
        _check_half_parent(parent, mspace_ndims(M), _symmetric_about_zero(parent))
        return new{M}(parent)
    end
end

@inline _check_half_parent(::Any, ::StaticInteger{0}, ::Any) = nothing
@inline _check_half_parent(μ, n::Integer, s) = n == 0 ? _check_half_parent(μ, static(0), s) : _throw_half_parent(μ, "scalar variates")
@inline _check_half_parent(μ, ::NoMSpaceElementSize, ::Any) = _throw_half_parent(μ, "a declared variate rank of zero")
@inline _check_half_parent(μ, ::StaticInteger{0}, ::False) = _throw_half_parent(μ, "symmetry about zero")
@noinline _throw_half_parent(μ, what) =
    throw(ArgumentError("Half requires measures with $what, but got a measure of type $(nameof(typeof(μ)))"))

# Whether a measure is symmetric about zero, `nothing` if not known:
_symmetric_about_zero(::Any) = nothing
_symmetric_about_zero(::Union{StdNormal,StdLogistic}) = static(true)
_symmetric_about_zero(::Union{StdExponential,StdUniform}) = static(false)
_symmetric_about_zero(μ::WeightedMeasure) = _symmetric_about_zero(μ.base)

@inline mspace_elsize(μ::Half) = mspace_elsize(μ.parent)
@inline mspace_flatsize(μ::Half) = mspace_flatsize(μ.parent)
@inline mspace_flatsize(::Type{<:Half{M}}) where {M} = mspace_flatsize(M)
@inline mspace_ndims(::Type{<:Half{M}}) where {M} = mspace_ndims(M)
@inline fixed_stream_size(::Type{<:Half{M}}) where {M} = fixed_stream_size(M)
@inline preferred_stdmeasure(::Type{<:Half}) = StdUniform

function Base.show(io::IO, μ::Half)
    print(io, "Half")
    show(io, μ.parent)
end

unhalf(μ::Half) = μ.parent

@inline function basemeasure(μ::Half)
    weightedmeasure(logtwo, basemeasure(unhalf(μ)))
end

@inline rand_impl(ctx::GenContext, μ::Half) = abs(rand_impl(ctx, unhalf(μ)))
@inline batched_rand_impl(ctx::GenContext, μ::Half, sz::SizeLike) = abs.(batched_rand_impl(ctx, unhalf(μ), sz))

function logdensityof_impl(μ::Half, x)
    ld = logdensityof(unhalf(μ), x) + logtwo
    return ifelse(x ≥ 0, ld, oftype(ld, -Inf))
end

logdensity_def(μ::Half, x) = logdensity_def(unhalf(μ), x)

@inline insupport(d::Half, x) = _insupport_and(x ≥ 0, insupport(unhalf(d), x))

testvalue(::Type{T}, ::Half) where {T} = one(T)

massof(μ::Half) = massof(unhalf(μ))

function smf(μ::Half, x)
    2 * smf(μ.parent, max(x, zero(x))) - 1
end

function invsmf(μ::Half, p)
    _nan_outside(StdUniform(), p, invsmf(μ.parent, _unit_interior((p + 1) / 2)))
end

@inline transport_to_std(::Type{StdUniform}, μ::Half, x) = _nan_outside(μ, x, smf(μ, x))
@inline transport_from_std(::Type{StdUniform}, μ::Half, p) = invsmf(μ, p)
