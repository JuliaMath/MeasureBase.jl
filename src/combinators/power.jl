import Base
using FillArrays: Fill
# """
# A power measure is a product of a measure with itself. The number of elements in
# the product determines the dimensionality of the resulting support.

# Note that power measures are only well-defined for integer powers.

# The nth power of a measure μ can be written μ^x.
# """
# PowerMeasure{M,N,D} = ProductMeasure{Fill{M,N,D}}

export PowerMeasure

struct PowerMeasure{M,A} <: AbstractProductMeasure
    parent::M
    axes::A
end

function Pretty.tile(μ::PowerMeasure)
    sz = length.(μ.axes)
    arg1 = Pretty.tile(μ.parent)
    arg2 = Pretty.tile(length(sz) == 1 ? only(sz) : sz)
    return Pretty.pair_layout(arg1, arg2; sep = " ^ ")
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::PowerMeasure{M}) where {T, M<:AbstractMeasure}
    map(CartesianIndices(d.axes)) do _
        rand(rng, T, d.parent)
    end
end

function Base.rand(rng::AbstractRNG, ::Type{T}, d::PowerMeasure) where {T}
    map(CartesianIndices(d.axes)) do _
        rand(rng, d.parent)
    end
end

@inline function powermeasure(x::T, sz::Tuple{Vararg{<:Any,N}}) where {T, N}
    a = axes(Fill{T, N}(x, sz))
    A = typeof(a)
    PowerMeasure{T,A}(x,a)
end

marginals(d::PowerMeasure) = Fill(d.parent, d.axes)

function Base.:^(μ::AbstractMeasure, dims::Tuple{Vararg{<:AbstractArray,N}}) where {N}
    powermeasure(μ, dims)
end

Base.:^(μ::AbstractMeasure, dims::Tuple) = powermeasure(μ, Base.OneTo.(dims))
Base.:^(μ::AbstractMeasure, n) = powermeasure(μ, (n,))

# Base.show(io::IO, d::PowerMeasure) = print(io, d.parent, " ^ ", size(d.xs))
# Base.show(io::IO, d::PowerMeasure{M,1}) where {M} = print(io, d.parent, " ^ ", length(d.xs))

# gentype(d::PowerMeasure{M,N}) where {M,N} = @inbounds Array{gentype(first(marginals(d))), N}

params(d::PowerMeasure) = params(first(marginals(d)))

# basemeasure(μ::PowerMeasure) = @inbounds basemeasure(first(μ.data))^size(μ.data)

@inline function basemeasure(d::PowerMeasure)
    basemeasure(d.parent) ^ d.axes
end

@inline function logdensity_def(d::PowerMeasure{M}, x) where {M}
    T = eltype(x)
    ℓ = 0.0
    # ℓ = zero(typeintersect(AbstractFloat,Core.Compiler.return_type(logdensity_def, Tuple{M,T})))
    parent = d.parent
    @inbounds for xj in x
        ℓ += logdensity_def(parent, xj)
    end
    ℓ
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