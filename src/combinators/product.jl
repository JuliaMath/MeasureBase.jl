using MappedArrays
using MappedArrays: ReadonlyMultiMappedArray
using Base: @propagate_inbounds
import Base
using FillArrays


"""
    abstract type AbstractProductMeasure

Abstact type for products of measures.

[`marginals(μ::AbstractProductMeasure)`](@ref) returns the collection of
measures that `μ` is the product of.
"""
abstract type AbstractProductMeasure <: AbstractMeasure end
export AbstractProductMeasure

function Pretty.tile(μ::AbstractProductMeasure)
    result = Pretty.literal("ProductMeasure(")
    result *= Pretty.tile(marginals(μ))
    result *= Pretty.literal(")")
end

massof(m::AbstractProductMeasure) = prod(massof, marginals(m))

function Base.:(==)(a::AbstractProductMeasure, b::AbstractProductMeasure)
    marginals(a) == marginals(b)
end
Base.length(μ::AbstractProductMeasure) = length(marginals(μ))
Base.size(μ::AbstractProductMeasure) = size(marginals(μ))


# TODO: Better `map` support in MappedArrays
_map(f, args...) = map(f, args...)
_map(f, x::MappedArrays.ReadonlyMappedArray) = mappedarray(fchain((x.f, f)), x.data)

function testvalue(::Type{T}, μ::AbstractProductMeasure) where {T}
    _map(m -> testvalue(T, m), marginals(μ))
end


function Base.rand(rng::AbstractRNG, ::Type{T}, d::AbstractProductMeasure) where {T}
    mar = marginals(d)
    _rand_product(rng, T, mar, eltype(mar))
end

function _rand_product(
    rng::AbstractRNG,
    ::Type{T},
    mar,
    ::Type{M},
) where {T,M<:AbstractMeasure}
    map(mar) do dⱼ
        rand(rng, T, dⱼ)
    end
end

function _rand_product(
    rng::AbstractRNG,
    ::Type{T},
    mar::ReadonlyMappedArray,
    ::Type{M},
) where {T,M<:AbstractMeasure}
    mappedarray(mar.data) do dⱼ
        rand(rng, T, mar.f(dⱼ))
    end |> collect
end

function _rand_product(rng::AbstractRNG, ::Type{T}, mar, ::Type{M}) where {T,M}
    map(mar) do dⱼ
        rand(rng, dⱼ)
    end
end

function _rand_product(
    rng::AbstractRNG,
    ::Type{T},
    mar::ReadonlyMappedArray,
    ::Type{M},
) where {T,M}
    mappedarray(mar.data) do dⱼ
        rand(rng, mar.f(dⱼ))
    end |> collect
end


@inline function logdensity_def(μ::AbstractProductMeasure, x)
    marginals_density_op(logdensity_def, marginals(μ), x)
end
@inline function unsafe_logdensityof(μ::AbstractProductMeasure, x)
    marginals_density_op(unsafe_logdensityof, marginals(μ), x)
end
@inline function logdensity_rel(μ::AbstractProductMeasure, ν::AbstractProductMeasure, x)
    marginals_density_op(logdensity_rel, marginals(μ), marginals(ν), x)
end

function _marginals_density_op(density_op::F, marginals_μ, x) where F
    mapreduce(density_op, +, marginals_μ, x)
end
@inline function _marginals_density_op(density_op::F, marginals_μ::Tuple, x::Tuple) where F
    # For tuples, `mapreduce` can have trouble with type inference
    sum(map(density_op, marginals_μ, x))
end
@inline function _marginals_density_op(density_op::F, marginals_μ::NameTuple{names}, x::NameTuple) where {F,names}
    nms = Val{names}()
    _marginals_density_op(density_op, marginals_μ, _reorder_nt(values(x), nms))
end

function _marginals_density_op(density_op::F, marginals_μ, marginals_ν, x) where F
    mapreduce(density_op, +, marginals_μ, marginals_ν, x)
end
@inline function _marginals_density_op(density_op::F, marginals_μ::Tuple, marginals_ν::Tuple, x::Tuple) where F
    # For tuples, `mapreduce` can have trouble with type inference
    sum(map(density_op, marginals_μ, marginals_ν, x))
end
@inline function _marginals_density_op(density_op::F, marginals_μ::NameTuple{names}, marginals_ν::NameTuple, x::NameTuple) where {F,names}
    nms = Val{names}()
    _marginals_density_op(density_op, marginals_μ, _reorder_nt(values(marginals_ν), nms), _reorder_nt(values(x), nms))
end


@inline basemeasure(μ::AbstractProductMeasure) =_marginals_basemeasure(marginals(μ))

_marginals_basemeasure(marginals_μ) = productmeasure(map(basemeasure, marginals_μ))

function _marginals_basemeasure(marginals_μ::Base.Generator{I,F}) where {I,F}
    T = Core.Compiler.return_type(marginals_μ.f, Tuple{eltype(marginals_μ.iter)})
    B = Core.Compiler.return_type(basemeasure, Tuple{T})
    _marginals_basemeasure_impl(μ, B, static(Base.issingletontype(B)))
end

function _marginals_basemeasure(marginals_μ::AbstractMappedArray{T}) where {T}
    B = Core.Compiler.return_type(basemeasure, Tuple{T})
    _marginals_basemeasure_impl(marginals_μ, B, static(Base.issingletontype(B)))
end

function _marginals_basemeasure_impl(marginals_μ, ::Type{B}, ::True) where {B}
    instance(B)^axes(marginals_μ)
end

function _marginals_basemeasure_impl(marginals_μ::AbstractMappedArray{T}, ::Type{B}, ::False) where {T,B}
    productmeasure(mappedarray(basemeasure, marginals_μ))
end

function _marginals_basemeasure_impl(marginals_μ::Base.Generator{I,F}, ::Type{B}, ::False) where {I,F,B}
    productmeasure(Base.Generator(basekernel(marginals_μ.f), marginals_μ.iter))
end



"""
    struct MeasureBase.ProductMeasure{M} <: AbstractProductMeasure

Represents a products of measures.

´ProductMeasure` wraps a collection of measures, this collection then
becomes the collection of the marginal measures of the `ProductMeasure`.

User code should not instantiate `ProductMeasure` directly, but should call
[`productmeasure`](@ref) instead.
"""
struct ProductMeasure{M} <: AbstractProductMeasure
    marginals::M
end

function Pretty.tile(d::ProductMeasure{T}) where {T<:Tuple}
    Pretty.list_layout(Pretty.tile.([marginals(d)...]), sep = " ⊗ ")
end

marginals(μ::ProductMeasure) = μ.marginals

proxy(μ::ProductMeasure{<:Fill}) = powermeasure(_fill_value(marginals(μ)), _fill_axes(marginals(μ)))


###############################################################################
# I <: Base.Generator

export rand!
using Random: rand!, GLOBAL_RNG, AbstractRNG

@propagate_inbounds function Random.rand!(
    rng::AbstractRNG,
    d::ProductMeasure,
    x::AbstractArray,
)
    # TODO: Generalize this
    T = Float64
    for (j, m) in zip(eachindex(x), marginals(d))
        @inbounds x[j] = rand(rng, T, m)
    end
    return x
end

export rand!
using Random: rand!, GLOBAL_RNG

function _rand(rng::AbstractRNG, ::Type{T}, d::ProductMeasure, mar::AbstractArray) where {T}
    elT = typeof(rand(rng, T, first(mar)))

    sz = size(mar)
    x = Array{elT,length(sz)}(undef, sz)
    rand!(rng, d, x)
end

@inline function insupport(d::AbstractProductMeasure, x::AbstractArray)
    mar = marginals(d)
    # We might get lucky and know statically that everything is inbounds
    T = Core.Compiler.return_type(insupport, Tuple{eltype(mar),eltype(x)})
    T <: True || all(zip(x, mar)) do (xj, mj)
        insupport(mj, xj) == true
    end
end

@inline function insupport(d::AbstractProductMeasure, x)
    for (mj, xj) in zip(marginals(d), x)
        dynamic(insupport(mj, xj)) || return false
    end
    return true
end

getdof(d::AbstractProductMeasure) = mapreduce(getdof, +, marginals(d))

function checked_arg(μ::ProductMeasure{<:NTuple{N,Any}}, x::NTuple{N,Any}) where {N}
    map(checked_arg, marginals(μ), x)
end

function checked_arg(
    μ::ProductMeasure{<:NamedTuple{names}},
    x::NamedTuple{names},
) where {names}
    NamedTuple{names}(map(checked_arg, values(marginals(μ)), values(x)))
end



