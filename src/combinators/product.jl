export ProductMeasure

using MappedArrays
using MappedArrays: ReadonlyMultiMappedArray
using Base: @propagate_inbounds
import Base
using FillArrays

export AbstractProductMeasure

abstract type AbstractProductMeasure <: AbstractMeasure end

function Pretty.tile(μ::AbstractProductMeasure)
    result = Pretty.literal("ProductMeasure(")
    result *= Pretty.tile(marginals(μ))
    result *= Pretty.literal(")")
end

massof(m::AbstractProductMeasure) = prod(massof, marginals(m))

export marginals

function Base.:(==)(a::AbstractProductMeasure, b::AbstractProductMeasure)
    marginals(a) == marginals(b)
end
Base.length(μ::AbstractProductMeasure) = length(marginals(μ))
Base.size(μ::AbstractProductMeasure) = size(marginals(μ))

basemeasure(d::AbstractProductMeasure) = productmeasure(map(basemeasure, marginals(d)))

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

@inline function logdensity_def(d::AbstractProductMeasure, x)
    mapreduce(logdensity_def, +, marginals(d), x)
end

struct ProductMeasure{M} <: AbstractProductMeasure
    marginals::M
end

@inline function logdensity_rel(μ::ProductMeasure, ν::ProductMeasure, x)
    mapreduce(logdensity_rel, +, marginals(μ), marginals(ν), x)
end

function Pretty.tile(d::ProductMeasure{T}) where {T<:Tuple}
    Pretty.list_layout(Pretty.tile.([marginals(d)...]), sep = " ⊗ ")
end

# For tuples, `mapreduce` has trouble with type inference
@inline function logdensity_def(d::ProductMeasure{T}, x) where {T<:Tuple}
    ℓs = map(logdensity_def, marginals(d), x)
    sum(ℓs)
end

@generated function logdensity_def(d::ProductMeasure{NamedTuple{N,T}}, x) where {N,T}
    k1 = QuoteNode(first(N))
    q = quote
        m = marginals(d)
        ℓ = logdensity_def(getproperty(m, $k1), getproperty(x, $k1))
    end
    for k in Base.tail(N)
        k = QuoteNode(k)
        qk = :(ℓ += logdensity_def(getproperty(m, $k), getproperty(x, $k)))
        push!(q.args, qk)
    end

    return q
end

# @generated function basemeasure(d::ProductMeasure{NamedTuple{N,T}}, x) where {N,T}
#     q = quote
#         m = marginals(d)
#     end
#     for k in N
#         qk = QuoteNode(k)
#         push!(q.args, :($k = basemeasure(getproperty(m, $qk))))
#     end

#     vals = map(x -> Expr(:(=), x,x), N)
#     push!(q.args, Expr(:tuple, vals...))
#     return q
# end

function basemeasure(μ::ProductMeasure{Base.Generator{I,F}}) where {I,F}
    mar = marginals(μ)
    T = Core.Compiler.return_type(mar.f, Tuple{eltype(mar.iter)})
    B = Core.Compiler.return_type(basemeasure, Tuple{T})
    _basemeasure(μ, B, static(Base.issingletontype(B)))
end

function basemeasure(μ::ProductMeasure{A}) where {T,A<:AbstractMappedArray{T}}
    B = Core.Compiler.return_type(basemeasure, Tuple{T})
    _basemeasure(μ, B, static(Base.issingletontype(B)))
end

function _basemeasure(μ::ProductMeasure, ::Type{B}, ::True) where {B}
    return instance(B)^axes(marginals(μ))
end

function _basemeasure(
    μ::ProductMeasure{A},
    ::Type{B},
    ::False,
) where {T,A<:AbstractMappedArray{T},B}
    mar = marginals(μ)
    productmeasure(mappedarray(basemeasure, mar))
end

function _basemeasure(
    μ::ProductMeasure{Base.Generator{I,F}},
    ::Type{B},
    ::False,
) where {I,F,B}
    mar = marginals(μ)
    productmeasure(Base.Generator(basekernel(mar.f), mar.iter))
end

marginals(μ::ProductMeasure) = μ.marginals

# TODO: Better `map` support in MappedArrays
_map(f, args...) = map(f, args...)
_map(f, x::MappedArrays.ReadonlyMappedArray) = mappedarray(fchain((x.f, f)), x.data)

function testvalue(::Type{T}, d::AbstractProductMeasure) where {T}
    _map(m -> testvalue(T, m), marginals(d))
end

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


function transport_to(ν::Pro)
end


function _marginal_transport_def(marginals_ν::NamedTuple{names}, marginals_μ::NamedTuple, x) where names
    # ToDo - Improvement: Match names as far as possible, even if in different order, and transport between
    # the rest in the order given.
    NamedTuple{names}(transport_to.(values(marginals_ν), values(marginals_μ), x))
end

@inline function _marginal_transport_def(marginals_ν, marginals_μ, x) 
    marginal_transport_non_ntnt(marginals_ν, marginals_μ, x)
end



function _marginal_transport_def(marginals_ν::AbstractVector{<:AbstractMeasure}, marginals_μ::AbstractVector{<:AbstractMeasure}, x)
    @assert x isa AbstractVector  # Sanity check, should not fail
    transport_to.(marginals_ν, marginals_μ, x)
end

function _marginal_transport_def(marginals_ν::Tuple{Vararg{AbstractMeasure,N}}, marginals_μ::Tuple{Vararg{AbstractMeasure,N}}, x) where N
    @assert x isa Tuple{Vararg{AbstractMeasure,N}}  # Sanity check, should not fail
    transport_to.(marginals_ν, marginals_μ, x)
end

function _marginal_transport_def(marginals_ν::NamedTuple{names}, marginals_μ::Tuple{Vararg{AbstractMeasure,N}}, x) where {names,N}
    _marginal_transport_def(marginals_ν, NamedTuple{names}(marginals_μ), x)
end

function _marginal_transport_def(marginals_ν::Tuple{Vararg{AbstractMeasure,N}}, marginals_μ::NamedTuple{names}, x) where {names,N}
    _marginal_transport_def(marginals_ν, values(marginals), x)
end

function _marginal_transport_def(marginals_ν::AbstractVector{<:AbstractMeasure}, marginals_μ::Tuple{Vararg{AbstractMeasure,N}}, x) where N
    _marginal_transport_def(_as_tuple(marginals_ν, Val(N)), marginals_μ, x)
end

function _marginal_transport_def(marginals_ν::Tuple{Vararg{AbstractMeasure,N}}, marginals_μ::AbstractVector{<:AbstractMeasure}, x) where N
    _marginal_transport_def(marginals_ν, _as_tuple(marginals_μ, Val(N)), x)
end



# Transport for products


#!!!!!!!!!!!!!!!!!!!!!! TODO:

# Helpers for product transforms and similar:

struct _TransportToStd{NU<:StdMeasure} <: Function end
_TransportToStd{NU}(μ, x) where {NU} = transport_to(NU()^getdof(μ), μ)(x)

struct _TransportFromStd{MU<:StdMeasure} <: Function end
_TransportFromStd{MU}(ν, x) where {MU} = transport_to(ν, MU()^getdof(ν))(x)

function _tuple_transport_def(
    ν::PowerMeasure{NU},
    μs::Tuple,
    xs::Tuple,
) where {NU<:StdMeasure}
    reshape(vcat(map(_TransportToStd{NU}, μs, xs)...), ν.axes)
end

function transport_def(
    ν::PowerMeasure{NU},
    μ::ProductMeasure{<:Tuple},
    x,
) where {NU<:StdMeasure}
    _tuple_transport_def(ν, marginals(μ), x)
end

function transport_def(
    ν::PowerMeasure{NU},
    μ::ProductMeasure{<:NamedTuple{names}},
    x,
) where {NU<:StdMeasure,names}
    _tuple_transport_def(ν, values(marginals(μ)), values(x))
end

@inline _offset_cumsum(s, x, y, rest...) = (s, _offset_cumsum(s + x, y, rest...)...)
@inline _offset_cumsum(s, x) = (s,)
@inline _offset_cumsum(s) = ()

function _stdvar_viewranges(μs::Tuple, startidx::IntegerLike)
    N = map(getdof, μs)
    offs = _offset_cumsum(startidx, N...)
    map((o, n) -> o:o+n-1, offs, N)
end

function _tuple_transport_def(
    νs::Tuple,
    μ::PowerMeasure{MU},
    x::AbstractArray{<:Real},
) where {MU<:StdMeasure}
    vrs = _stdvar_viewranges(νs, firstindex(x))
    xs = map(r -> view(x, r), vrs)
    map(_TransportFromStd{MU}, νs, xs)
end

function transport_def(
    ν::ProductMeasure{<:Tuple},
    μ::PowerMeasure{MU},
    x,
) where {MU<:StdMeasure}
    _tuple_transport_def(marginals(ν), μ, x)
end

function transport_def(
    ν::ProductMeasure{<:NamedTuple{names}},
    μ::PowerMeasure{MU},
    x,
) where {MU<:StdMeasure,names}
    NamedTuple{names}(_tuple_transport_def(values(marginals(ν)), μ, x))
end
