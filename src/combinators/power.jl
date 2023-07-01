import Base

"""
    marginals(μ::AbstractMeasure)

Returns the marginals measures of `μ` as a collection of measures.

The kind of marginalization implied by `marginals` depends on the
type of `μ`.

`μ` may be a power of a measure or a product of measures, but other
types of measures may support `marginals` as well.
"""
function marginals end
export marginals


export PowerMeasure

"""
    struct PowerMeasure{M,...} <: AbstractProductMeasure

A power measure is a product of a measure with itself. The number of elements in
the product determines the dimensionality of the resulting support.

Note that power measures are only well-defined for integer powers.

The nth power of a measure μ can be written μ^n.
"""
struct PowerMeasure{M,A} <: AbstractProductMeasure
    parent::M
    axes::A
end

maybestatic_length(μ::PowerMeasure) = prod(maybestatic_size(μ))
maybestatic_size(μ::PowerMeasure) = map(maybestatic_length, μ.axes)

function Pretty.tile(μ::PowerMeasure)
    sz = length.(μ.axes)
    arg1 = Pretty.tile(μ.parent)
    arg2 = Pretty.tile(length(sz) == 1 ? only(sz) : sz)
    return Pretty.pair_layout(arg1, arg2; sep = " ^ ")
end

# ToDo: Make rand return static arrays for statically-sized power measures.

function _cartidxs(axs::Tuple{Vararg{AbstractUnitRange,N}}) where {N}
    CartesianIndices(map(_dynamic, axs))
end

function Base.rand(
    rng::AbstractRNG,
    ::Type{T},
    d::PowerMeasure{M},
) where {T,M<:AbstractMeasure}
    map(_cartidxs(d.axes)) do _
        rand(rng, T, d.parent)
    end
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::PowerMeasure) where {T}
    map(_cartidxs(d.axes)) do _
        rand(rng, d.parent)
    end
end

@inline _pm_axes(sz::Tuple{Vararg{IntegerLike,N}}) where {N} = map(one_to, sz)
@inline _pm_axes(axs::Tuple{Vararg{AbstractUnitRange,N}}) where {N} = axs

@inline function powermeasure(x::T, sz::Tuple{Vararg{Any,N}}) where {T,N}
    PowerMeasure(x, _pm_axes(sz))
end

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

@inline function logdensity_def(d::PowerMeasure{M}, x) where {M}
    parent = d.parent
    sum(x) do xj
        logdensity_def(parent, xj)
    end
end

@inline function logdensity_def(d::PowerMeasure{M,Tuple{Static.SOneTo{N}}}, x) where {M,N}
    parent = d.parent
    sum(1:N) do j
        @inbounds logdensity_def(parent, x[j])
    end
end

@inline function logdensity_def(
    d::PowerMeasure{M,NTuple{N,Static.SOneTo{0}}},
    x,
) where {M,N}
    static(0.0)
end

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

@inline getdof(μ::PowerMeasure) = getdof(μ.parent) * prod(map(length, μ.axes))

@inline function getdof(::PowerMeasure{<:Any,NTuple{N,Static.SOneTo{0}}}) where {N}
    static(0)
end

@propagate_inbounds function checked_arg(μ::PowerMeasure, x::AbstractArray{<:Any})
    @boundscheck begin
        sz_μ = map(length, μ.axes)
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

# To avoid ambiguities
function logdensity_def(
    ::PowerMeasure{P,Tuple{Vararg{Static.SOneTo{0},N}}},
    x,
) where {P<:PrimitiveMeasure,N}
    static(0.0)
end


# For transport, always pull back to one-dimensional PowerMeasure first:

transport_origin(μ::PowerMeasure{<:Any,N}) where N = ν.parent^product(map(length, μ.axes))

function from_origin(μ::_PowerStdMeasure{<:Any,N}, x_origin) where N
    # Sanity check, should never fail:
    @assert x_origin isa AbstractVector
    return reshape(x_origin, map(length, μ.axes)...)
end


# One-dimensional PowerMeasure has an origin iff it's parent has an origin:

transport_origin(μ::PowerMeasure{<:AbstractMeasure,1}) = _origin_pwr(::typeof(μ), transport_origin(μ.parent), μ.axes)
_pwr_origin(::Type{MU}, parent_origin, axes) = parent_origin^axes
_pwr_origin(::Type{MU}, ::NoTransportOrigin, axes) = NoTransportOrigin{MU}

function from_origin(μ::PowerMeasure{<:AbstractMeasure,1}, x_origin)
    # Sanity check, should never fail:
    @assert x_origin isa AbstractVector
    from_origin.(Ref(μ.parent), x_origin)
end

