import Base


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

maybestatic_length(μ::PowerMeasure) = prod(maybestatic_size(μ))
maybestatic_size(μ::PowerMeasure) = map(maybestatic_length, μ.axes)


"""
    MeasureBase.pwr_base(μ::PowerMeasure)

Returns `ν` for `μ = ν^axs`
"""
pwr_base(μ::PowerMeasure) = μ.parent


"""
    MeasureBase.pwr_axes(μ::PowerMeasure)

Returns `axs` for `μ = ν^axs`, `axs` being a tuple of integer ranges.
"""
pwr_axes(μ::PowerMeasure) = μ.axes


"""
    MeasureBase.pwr_size(μ::PowerMeasure)

Returns `sz` for `μ = ν^sz`, `sz` being a tuple of integers.
"""
pwr_size(μ::PowerMeasure) = map(maybestatic_length, μ.axes)


function Pretty.tile(μ::PowerMeasure)
    sz = length.(μ.axes)
    arg1 = Pretty.tile(μ.parent)
    arg2 = Pretty.tile(length(sz) == 1 ? only(sz) : sz)
    return Pretty.pair_layout(arg1, arg2; sep = " ^ ")
end

# ToDo: Make rand and testvalue return static arrays for statically-sized power measures.

function _cartidxs(axs::Tuple{Vararg{AbstractUnitRange,N}}) where {N}
    CartesianIndices(map(_dynamic, axs))
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::PowerMeasure) where {T}
    map(_ -> rand(rng, T, d.parent), _cartidxs(d.axes))
end

function Base.rand(rng::AbstractRNG, d::PowerMeasure)
    map(_ -> rand(rng, d.parent), _cartidxs(d.axes))
end

function testvalue(::Type{T}, d::PowerMeasure) where {T}
    map(_ -> testvalue(T, d.parent), _cartidxs(d.axes))
end

function testvalue(d::PowerMeasure)
    map(_ -> testvalue(d.parent), _cartidxs(d.axes))
end


@inline _pm_axes(sz::Tuple{Vararg{IntegerLike,N}}) where {N} = map(one_to, sz)
@inline _pm_axes(axs::Tuple{Vararg{AbstractUnitRange,N}}) where {N} = axs

marginals(d::PowerMeasure) = fill_with(d.parent, d.axes)

function Base.:^(μ::AbstractMeasure, dims::Tuple{Vararg{AbstractArray,N}}) where {N}
    powermeasure(μ, dims)
end

Base.:^(μ::AbstractMeasure, dims::Tuple) = powermeasure(μ, one_to.(dims))
Base.:^(μ::AbstractMeasure, n) = powermeasure(μ, (n,))

# Base.show(io::IO, d::PowerMeasure) = print(io, d.parent, " ^ ", size(d.xs))
# Base.show(io::IO, d::PowerMeasure{M,1}) where {M} = print(io, d.parent, " ^ ", length(d.xs))

# gentype(d::PowerMeasure{M,N}) where {M,N} = @inbounds Array{gentype(first(marginals(d))), N}

params(d::PowerMeasure) = params(first(marginals(d)))

# basemeasure(μ::PowerMeasure) = @inbounds basemeasure(first(μ.data))^size(μ.data)

@inline function basemeasure(d::PowerMeasure)
    basemeasure(d.parent)^d.axes
end

@inline logdensity_def(d::PowerMeasure, x) = _pwr_logdensity_def(d.parent, x, prod(pwr_size(d)))

@inline _pwr_logdensity_def(::PowerMeasure, x, ::Integer, ::StaticInteger{0}) = static(false)

@inline function _pwr_logdensity_def(d::PowerMeasure, x, ::IntegerLike)
    parent = d.parent
    sum(x) do xj
        logdensity_def(parent, xj)
    end
end

# ToDo: Specialized version of _pwr_logdensity_def for statically-sized power measures

# ToDo: Re-enable this?
# _pwr_logdensity_def(::PowerMeasure{P}, x, ::IntegerLike) where {P<:PrimitiveMeasure} = static(0.0)


@inline function insupport(μ::PowerMeasure, x)
    p = μ.parent
    all(x) do xj
        # https://github.com/SciML/Static.jl/issues/36
        dynamic(insupport(p, xj))
    end
end

@inline function insupport(μ::PowerMeasure, x::AbstractArray)
    p = μ.parent
    all(x) do xj
        # https://github.com/SciML/Static.jl/issues/36
        dynamic(insupport(p, xj))
    end
end

@inline getdof(μ::PowerMeasure) = getdof(μ.parent) * prod(pwr_size(μ))
@inline fast_dof(μ::PowerMeasure) = fast_dof(μ.parent) * prod(pwr_size(μ))

@inline function getdof(::PowerMeasure{<:Any,NTuple{N,Static.SOneTo{0}}}) where {N}
    static(0)
end

@inline function fast_dof(::PowerMeasure{<:Any,NTuple{N,Static.SOneTo{0}}}) where {N}
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


"""
    MeasureBase.StdPowerMeasure{MU<:StdMeasure,N}

Represents and N-dimensional power of the standard measure `MU()`.
"""
const StdPowerMeasure{MU<:StdMeasure,N} = PowerMeasure{MU,<:NTuple{N,Base.OneTo}}
