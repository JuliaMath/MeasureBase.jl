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

# Powers of scalar-variate measures have array-valued variates of known size:
@inline mspace_elsize(μ::PowerMeasure) = _pwr_mspace_elsize(μ, mspace_elsize(pwr_base(μ)))
@inline _pwr_mspace_elsize(μ::PowerMeasure, ::Tuple{}) = pwr_size(μ)
@inline _pwr_mspace_elsize(μ::PowerMeasure, ::Any) = NoMSpaceElementSize{typeof(μ)}()

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

# Power structure is unwrapped into the power axes arguments of the batched
# density machinery (see density-batched.jl), which fuses evaluation over
# flat variate storage:

for head in [:logdensityof_impl, :logdensity_def]
    @eval @inline function $head(d::PowerMeasure, x)
        _powered_ld($head, pwr_base(d), x, pwr_axes(d))
    end

    @eval @inline function $head(
        ::PowerMeasure{<:Any,<:Tuple{Vararg{StaticOneToLike{0}}}},
        x,
    )
        static(0.0)
    end
end

@inline function insupport(μ::PowerMeasure, x)
    p = μ.parent
    all(x) do xj
        # https://github.com/SciML/Static.jl/issues/36
        dynamic(insupport(p, xj))
    end
end

_all(A) = all(A)
_all(::AbstractArray{NoFastInsupport{T}}) where {T} = NoFastInsupport{T}()

@inline function insupport(μ::PowerMeasure, x::AbstractArray)
    p = μ.parent
    insupp = broadcast(x) do xj
        # https://github.com/SciML/Static.jl/issues/36
        dynamic(insupport(p, xj))
    end
    _all(insupp)
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

@propagate_inbounds function checked_arg(μ::PowerMeasure, x::AbstractArray{<:Any})
    @boundscheck begin
        sz_μ = pwr_size(μ)
        sz_x = size(x)
        if sz_μ != sz_x
            throw(ArgumentError("Size of variate doesn't match size of power measure"))
        end
    end
    return x
end

function checked_arg(μ::PowerMeasure, x::Any)
    throw(ArgumentError("Size of variate doesn't match size of power measure"))
end

massof(m::PowerMeasure) = massof(m.parent)^prod(m.axes)

logdensity_def(::PowerMeasure{P}, x) where {P<:PrimitiveMeasure} = static(0.0)

# Disambiguation with the static-zero-size power density method:
function logdensity_def(
    ::PowerMeasure{P,<:Tuple{Vararg{StaticOneToLike{0}}}},
    ::Any,
) where {P<:PrimitiveMeasure}
    static(0.0)
end

