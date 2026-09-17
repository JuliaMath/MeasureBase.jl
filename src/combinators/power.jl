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

# ToDo: Make rand return static arrays for statically-sized power measures.

function _cartidxs(axs::Tuple{Vararg{AbstractUnitRange,N}}) where {N}
    CartesianIndices(map(asnonstatic, axs))
end

function Base.rand(
    rng::AbstractRNG,
    ::Type{T},
    d::PowerMeasure{M},
) where {T,M<:AbstractMeasure}
    axs, base_d = pwr_axes(d), pwr_base(d)
    map(_cartidxs(axs)) do _
        rand(rng, T, base_d)
    end
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::PowerMeasure) where {T}
    axs, base_d = pwr_axes(d), pwr_base(d)
    map(_cartidxs(axs)) do _
        rand(rng, base_d)
    end
end

marginals(d::PowerMeasure) = maybestatic_fill(d.parent, d.axes)

@inline mspace_elsize(μ::PowerMeasure) = pwr_size(μ)
@inline mspace_flatsize(μ::PowerMeasure) = _cat_sizes(mspace_flatsize(pwr_base(μ)), pwr_size(μ))

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

# Densities of powers are evaluated by the batched density machinery over
# the flat variate storage (see density-batched.jl):

@inline logdensityof_impl(μ::PowerMeasure, x) = _powered_ld(logdensityof_impl, μ, x)
@inline logdensity_def(μ::PowerMeasure, x) = _powered_ld(logdensity_def, μ, x)
@inline batched_logdensityof_impl(μ::PowerMeasure, A::AbstractArray) = _batched_ld(logdensityof_impl, μ, A)

# Support checks of powers run over the flat variate storage where the base
# measure has scalar variates, elementwise otherwise:
@inline function insupport(μ::PowerMeasure, x::AbstractArray)
    _powered_insupport(μ, x, _flat_storage(x), mspace_flatsize(μ))
end

@inline function _powered_insupport(μ::PowerMeasure, x, x_flat::AbstractArray, ::SizeLike)
    ν, _ = _pwr_unwrap(μ)
    _powered_insupport_flat(ν, x_flat, mspace_flatsize(ν))
end
@inline function _powered_insupport_flat(ν, x_flat::AbstractArray, ::Tuple{})
    _all_insupport(broadcast(_insupport_bool ∘ Base.Fix1(insupport, ν), x_flat))
end
@inline _powered_insupport_flat(ν, x_flat::AbstractArray, ::Any) = _powered_insupport_elementwise(ν, x_flat)
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
    @boundscheck begin
        sz_x = maybestatic_size(x)
        if sz_x != pwr_size(μ) && !_matches_flatsize(sz_x, mspace_flatsize(μ))
            _throw_size_mismatch()
        end
    end
    return x
end

@inline _matches_flatsize(sz_x, sz_flat::SizeLike) = Tuple(sz_x) == Tuple(sz_flat)
@inline _matches_flatsize(sz_x, ::NoMSpaceElementSize) = false

checked_arg(μ::PowerMeasure, x::Any) = _throw_size_mismatch()

massof(m::PowerMeasure) = massof(m.parent)^prod(m.axes)


