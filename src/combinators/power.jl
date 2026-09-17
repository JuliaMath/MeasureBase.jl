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


# Transport: the standard variate of a power is the flat vector of the
# standard variates of its base measure, in the order of the flat variate
# storage.

function transport_to_std(::Type{S}, μ::PowerMeasure, x::AbstractArray) where {S<:StdMeasure}
    _pwr_to_std(S, μ, x, _flat_storage(x), mspace_flatsize(μ))
end

# Flat storage of known flat size: transport the variates of the innermost
# base measure over the flat storage.
function _pwr_to_std(::Type{S}, μ::PowerMeasure, x::AbstractArray, x_flat::AbstractArray, sz_flat::SizeLike) where {S}
    ν, _ = _pwr_unwrap(μ)
    _check_flatsize(x_flat, sz_flat)
    _pwr_to_std_flat(S, ν, x_flat, mspace_flatsize(ν))
end

@inline function _pwr_to_std_flat(::Type{S}, ν, x_flat::AbstractArray, ::Tuple{}) where {S}
    _flat_std_of(broadcast(Base.Fix1(_ToStd{S}(), ν), x_flat))
end

@inline function _pwr_to_std_flat(::Type{S}, ν, x_flat::AbstractArray, sz::SizeLike) where {S}
    _flat_std_of(map(Base.Fix1(_ToStd{S}(), ν), sliced(x_flat, Val(length(sz)))))
end

# Otherwise transport the variates of the base measure one by one:
function _pwr_to_std(::Type{S}, μ::PowerMeasure, x::AbstractArray, ::Any, ::Any) where {S}
    _flat_std_of(map(Base.Fix1(_ToStd{S}(), pwr_base(μ)), x))
end

function transport_from_std(::Type{S}, μ::PowerMeasure, z::AbstractVector) where {S<:StdMeasure}
    _check_stdlength(z, fast_dof(μ))
    _pwr_from_std(S, μ, z, mspace_flatsize(μ))
end

function _pwr_from_std(::Type{S}, μ::PowerMeasure, z::AbstractVector, sz_flat::SizeLike) where {S}
    ν, _ = _pwr_unwrap(μ)
    _pwr_variate(μ, _pwr_from_std_flat(S, ν, z, sz_flat, mspace_flatsize(ν)))
end

# Base measures of unknown variate size are transported one by one:
function _pwr_from_std(::Type{S}, μ::PowerMeasure, z::AbstractVector, ::NoMSpaceElementSize) where {S}
    ys, z_rest = _marginals_from_std_with_rest(S, marginals(μ), z)
    if !isempty(z_rest)
        throw(ArgumentError("Length of standard variate doesn't match degrees of freedom of power measure"))
    end
    return ys
end

@inline function _check_stdlength(z::AbstractVector, n::IntegerLike)
    if maybestatic_length(z) != n
        throw(ArgumentError("Length of standard variate doesn't match degrees of freedom of measure"))
    end
    return nothing
end
@inline _check_stdlength(::AbstractVector, ::AbstractNoDOF) = nothing

@inline function _pwr_from_std_flat(::Type{S}, ν, z::AbstractVector, sz_flat, ::Tuple{}) where {S}
    maybestatic_reshape(broadcast(Base.Fix1(_FromStd{S}(), ν), z), sz_flat)
end

@inline function _pwr_from_std_flat(::Type{S}, ν, z::AbstractVector, sz_flat, sz_ν::SizeLike) where {S}
    n_variates = size2length(sz_flat) ÷ size2length(sz_ν)
    maybestatic_reshape(stacked(_pwr_from_std_chunks(S, ν, z, n_variates, fast_dof(ν))), sz_flat)
end

function _pwr_from_std_chunks(::Type{S}, ν, z::AbstractVector, n_variates, dof_ν::IntegerLike) where {S}
    chunks = sliced(maybestatic_reshape(z, (dof_ν, n_variates)), Val(1))
    map(Base.Fix1(_FromStd{S}(), ν), chunks)
end

function _pwr_from_std_chunks(::Type{S}, ν, z::AbstractVector, n_variates, ::AbstractNoDOF) where {S}
    ys, z_rest = _marginals_from_std_with_rest(S, FillArrays.Fill(ν, n_variates), z)
    if !isempty(z_rest)
        throw(ArgumentError("Length of standard variate doesn't match degrees of freedom of power measure"))
    end
    return ys
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

# The nested variate layout of a power over its flat storage:
@inline _pwr_variate(μ::PowerMeasure, A::AbstractArray) = _pwr_nest(pwr_base(μ), _pwr_variate(pwr_base(μ), A))
@inline _pwr_variate(ν, A::AbstractArray) = _nest_leaf(A, mspace_flatsize(ν))
@inline _nest_leaf(A::AbstractArray, ::Tuple{}) = A
@inline _nest_leaf(A::AbstractArray, ::NoMSpaceElementSize) = A
@inline _nest_leaf(A::AbstractArray{<:Any,N}, sz::SizeLike) where {N} = _nest_leaf(A, Val(length(sz)), Val(N))
@inline _nest_leaf(A::AbstractArray, ::Val{N}, ::Val{N}) where {N} = A
@inline _nest_leaf(A::AbstractArray, ::Val{M}, ::Val) where {M} = sliced(A, Val(M))
@inline _pwr_nest(ν::PowerMeasure, B::AbstractArray) = sliced(B, Val(length(pwr_axes(ν))))
@inline _pwr_nest(ν, B::AbstractArray) = B
