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

rand_impl(ctx::GenContext, d::AbstractProductMeasure) = map(Base.Fix1(_marginal_rand, ctx), marginals(d))

@inline _marginal_rand(ctx::GenContext, m::AbstractMeasure) = rand_impl(ctx, m)
@inline _marginal_rand(ctx::GenContext, d) = rand(get_rng(ctx), d)

for (head, func) in [(:logdensityof_impl, :logdensityof), (:logdensity_def, :logdensity_def)]
    @eval @inline function $head(d::AbstractProductMeasure, x)
        _check_marginal_count(marginals(d), x)
        mapreduce($func, +, marginals(d), x)
    end
end

# Variates of products are collections of marginal variates, with the same
# structure as the marginals:
@inline function _check_marginal_count(mar::AbstractArray, x::AbstractArray)
    size(mar) == size(x) || _throw_marginal_mismatch()
    return nothing
end
@inline _check_marginal_count(::AbstractArray, x) = _throw_marginal_mismatch()
# Tuple products also take vector variates (e.g. from converted product
# distributions):
@inline function _check_marginal_count(mar::Tuple, x::Union{Tuple,AbstractVector})
    length(mar) == length(x) || _throw_marginal_mismatch()
    return nothing
end
@inline _check_marginal_count(::Tuple, x) = _throw_marginal_mismatch()
@inline _check_marginal_count(mar, x) = nothing

@noinline _throw_marginal_mismatch() =
    throw(ArgumentError("Variate doesn't match the structure of the marginals of a product measure"))

struct ProductMeasure{M} <: AbstractProductMeasure
    marginals::M
end

proxy(μ::ProductMeasure{<:FillArrays.Fill}) =
    powermeasure(_fill_value(marginals(μ)), _fill_axes(marginals(μ)))

# Relative densities between products evaluate marginal-wise. Support
# checks happen at the logdensity_rel level for the whole products, so the
# unsafe marginal evaluation suffices here:
@inline function logdensity_rel_def(μ::ProductMeasure, ν::ProductMeasure, x)
    mapreduce(unsafe_logdensity_rel, +, marginals(μ), marginals(ν), x)
end

# For tuples, `mapreduce` has trouble with type inference:
@inline function logdensity_rel_def(
    μ::ProductMeasure{<:Tuple},
    ν::ProductMeasure{<:Tuple},
    x,
)
    sum(map(unsafe_logdensity_rel, marginals(μ), marginals(ν), x))
end

_mspace_names(μ::ProductMeasure{<:NamedTuple{names}}) where {names} = names

function Pretty.tile(d::ProductMeasure{T}) where {T<:Tuple}
    Pretty.list_layout(Pretty.tile.([marginals(d)...]), sep = " ⊗ ")
end

@eval @generated function _product_gen_impl(
    ::Val{func},
    d::ProductMeasure{NamedTuple{N,T}},
    x,
) where {func,N,T}
    k1 = QuoteNode(first(N))
    q = quote
        m = marginals(d)
        ℓ = $func(getproperty(m, $k1), getproperty(x, $k1))
    end
    for k in Base.tail(N)
        k = QuoteNode(k)
        qk = :(ℓ += $func(getproperty(m, $k), getproperty(x, $k)))
        push!(q.args, qk)
    end

    return q
end

for (head, func) in [(:logdensityof_impl, :logdensityof), (:logdensity_def, :logdensity_def)]
    # For tuples, `mapreduce` has trouble with type inference
    @eval @inline function $head(d::ProductMeasure{T}, x) where {T<:Tuple}
        _check_marginal_count(marginals(d), x)
        ℓs = map($func, marginals(d), x)
        sum(ℓs)
    end

    @eval function $head(d::ProductMeasure{NamedTuple{N,T}}, x) where {N,T}
        _product_gen_impl(Val($func), d, x)
    end
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

"""
    MeasureBase.basekernel(f)

For a function `f` that returns a measure, return the function that returns
the base measure instead, satisfying `basekernel(f)(p) == basemeasure(f(p))`.
"""
function basekernel end

basekernel(f) = basemeasure ∘ f
basekernel(f::Returns) = Returns(basemeasure(f.value))

function _basemeasure(
    μ::ProductMeasure{Base.Generator{I,F}},
    ::Type{B},
    ::False,
) where {I,F,B}
    mar = marginals(μ)
    productmeasure(Base.Generator(basekernel(mar.f), mar.iter))
end

marginals(μ::ProductMeasure) = μ.marginals

@inline mspace_elsize(μ::ProductMeasure{<:AbstractArray}) = maybestatic_size(marginals(μ))

@inline function mspace_flatsize(μ::ProductMeasure{<:AbstractArray{M}}) where {M}
    _cat_sizes(mspace_flatsize(M), maybestatic_size(marginals(μ)))
end

# Batched densities over flat storage `(marginal flat dims..., product
# dims..., batch dims...)`. Marginals with scalar variates align with the
# leading dimensions of the batch, so one broadcast evaluates all marginal
# densities. Marginals with array variates are evaluated one by one over
# their slices of the batch.
@inline function batched_logdensityof_impl(μ::ProductMeasure{<:AbstractArray{M,N}}, A::AbstractArray) where {M,N}
    _product_batched_ld(μ, A, mspace_flatsize(M), Val(N), Val(isconcretetype(M)))
end

@inline function _product_batched_ld(μ::ProductMeasure, A::AbstractArray, ::Tuple{}, ::Val{N}, ::Val{true}) where {N}
    ld = Broadcast.instantiate(Broadcast.broadcasted(dynamic ∘ logdensityof_impl, marginals(μ), A))
    _sum_leading_dims(ld, static(N))
end

@inline function _product_batched_ld(μ::ProductMeasure, A::AbstractArray, sz_m::SizeLike, ::Val{N}, ::Val{true}) where {N}
    _marginal_slices_ld(marginals(μ), A, Val(length(sz_m)), Val(ndims(A) - length(sz_m) - N))
end

@inline function _product_batched_ld(μ::ProductMeasure, A::AbstractArray, ::Any, ::Val, ::Val)
    _batched_ld_generic(logdensityof_impl, μ, A)
end

function _marginal_slices_ld(mar::AbstractArray{<:Any,N}, A::AbstractArray, ::Val{K}, ::Val{B}) where {N,K,B}
    if ndims(A) != K + N + B || ntuple(i -> size(A, K + i), Val(N)) != size(mar)
        _throw_size_mismatch()
    end
    lead = ntuple(_ -> Colon(), Val(K))
    trail = ntuple(_ -> Colon(), Val(B))
    ld(i) = _materialize(batched_logdensityof_impl(mar[i], view(A, lead..., Tuple(i)..., trail...)))
    init = _zero_logd(A, ntuple(i -> size(A, K + N + i), Val(B)))
    return mapreduce(ld, +, CartesianIndices(mar); init = init)
end

@inline _zero_logd(A::AbstractArray, ::Tuple{}) = zero(_logd_numtype(A))
@inline _zero_logd(A::AbstractArray, dims::Tuple) = fill!(similar(A, _logd_numtype(A), dims), 0)

# The point density of array products with array-variate marginals accepts
# the flat variate storage `(marginal flat dims..., product dims...)`:
@inline function logdensityof_impl(μ::ProductMeasure{<:AbstractArray{M}}, x::AbstractArray{<:Number}) where {M}
    _array_product_ld(μ, x, mspace_flatsize(M))
end
@inline _array_product_ld(μ::ProductMeasure, x::AbstractArray, ::Tuple{}) = _array_product_ld_nested(μ, x)
@inline _array_product_ld(μ::ProductMeasure, x::AbstractArray, ::NoMSpaceElementSize) = _array_product_ld_nested(μ, x)
@inline function _array_product_ld(μ::ProductMeasure, x::AbstractArray, sz_m::SizeLike)
    _marginal_slices_ld(marginals(μ), x, Val(length(sz_m)), Val(0))
end
@inline function _array_product_ld_nested(μ::ProductMeasure, x::AbstractArray)
    _check_marginal_count(marginals(μ), x)
    mapreduce(logdensityof, +, marginals(μ), x)
end

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

@inline function insupport(d::AbstractProductMeasure, x::AbstractArray)
    _all_insupport(broadcast(_insupport_bool ∘ insupport, marginals(d), x))
end

@inline function insupport(d::AbstractProductMeasure, x)
    mapreduce(insupport, _insupport_and, marginals(d), x)
end

@inline _all_insupport(A::AbstractArray{<:NoFastInsupport{T}}) where {T} = NoFastInsupport{T}()
@inline _all_insupport(A::AbstractArray) = all(A)

getdof(d::AbstractProductMeasure) = _sum_dofs(getdof, marginals(d))
fast_dof(d::AbstractProductMeasure) = _sum_dofs(fast_dof, marginals(d))

# Sums over static DOFs of tuples fold at compile time, arrays of marginals
# are summed dynamically (also on GPU arrays):
@inline _sum_dofs(f, mar) = sum(f, mar)
@inline _sum_dofs(f, mar::AbstractArray) = mapreduce(_dynamic_dof ∘ f, +, mar; init = 0)
@inline _dynamic_dof(n::IntegerLike) = dynamic(n)
@inline _dynamic_dof(nodof::AbstractNoDOF) = nodof

function checked_arg(μ::ProductMeasure{<:NTuple{N,Any}}, x::NTuple{N,Any}) where {N}
    map(checked_arg, marginals(μ), x)
end

# Variates of array products are arrays of marginal variates or, for
# marginals with array variates of known size, their flat storage:
@propagate_inbounds function checked_arg(μ::ProductMeasure{<:AbstractArray{M}}, x::AbstractArray) where {M}
    @boundscheck _check_product_arg(marginals(μ), x, mspace_flatsize(M))
    return x
end

@inline _check_product_arg(mar, x::AbstractArray, ::Tuple{}) = _check_marginal_count(mar, x)
@inline _check_product_arg(mar, x::AbstractArray{<:Number}, ::NoMSpaceElementSize) = _check_marginal_count(mar, x)
@inline function _check_product_arg(mar, x::AbstractArray, ::NoMSpaceElementSize)
    _check_marginal_count(mar, x)
    foreach(checked_arg, mar, x)
    return nothing
end
@inline function _check_product_arg(mar, x::AbstractArray, sz_m::SizeLike)
    if size(x) == size(mar)
        foreach(checked_arg, mar, x)
    elseif size(x) != (Tuple(sz_m)..., size(mar)...)
        _throw_marginal_mismatch()
    end
    return nothing
end


function checked_arg(
    μ::ProductMeasure{<:NamedTuple{names}},
    x::NamedTuple{names},
) where {names}
    NamedTuple{names}(map(checked_arg, values(marginals(μ)), values(x)))
end


# Transport marginal by marginal, the standard variates of the marginals
# are concatenated in order:

function transport_to_std(::Type{S}, μ::ProductMeasure{<:Tuple}, x::Tuple) where {S<:StdMeasure}
    _flatten_to_rv(map((m, xi) -> _as_stdstream(transport_to_std(S, m, xi)), marginals(μ), x))
end

function transport_to_std(::Type{S}, μ::ProductMeasure{<:NamedTuple{names}}, x::NamedTuple{names}) where {S<:StdMeasure,names}
    transport_to_std(S, productmeasure(values(marginals(μ))), values(x))
end

function transport_to_std(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, x::AbstractArray) where {S<:StdMeasure,M}
    _array_product_to_std(S, μ, x, Val(isconcretetype(M)))
end
function _array_product_to_std(::Type{S}, μ, x::AbstractArray, ::Val{true}) where {S}
    _flat_std_of(broadcast(_ToStd{S}(), marginals(μ), x))
end
# Marginals of mixed types may have standard variates of mixed shapes:
function _array_product_to_std(::Type{S}, μ, x::AbstractArray, ::Val{false}) where {S}
    zs = [_as_stdstream(transport_to_std(S, m, xi)) for (m, xi) in zip(marginals(μ), x)]
    isempty(zs) ? SVector{0,Bool}() : reduce(vcat, zs)
end

function transport_from_std_with_rest(::Type{S}, μ::ProductMeasure{<:Tuple}, z::AbstractVector) where {S<:StdMeasure}
    _marginals_from_std_with_rest(S, marginals(μ), z)
end

function transport_from_std_with_rest(::Type{S}, μ::ProductMeasure{<:NamedTuple{names}}, z::AbstractVector) where {S<:StdMeasure,names}
    ys, z_rest = _marginals_from_std_with_rest(S, values(marginals(μ)), z)
    return NamedTuple{names}(ys), z_rest
end

function transport_from_std_with_rest(::Type{S}, μ::ProductMeasure{<:AbstractArray}, z::AbstractVector) where {S<:StdMeasure}
    _array_product_from_std_with_rest(S, μ, z, fast_dof(μ))
end
@inline function _array_product_from_std_with_rest(::Type{S}, μ, z, n::IntegerLike) where {S}
    _from_std_with_rest_bydof(S, μ, z, n)
end
function _array_product_from_std_with_rest(::Type{S}, μ, z, ::AbstractNoDOF) where {S}
    _marginals_from_std_with_rest(S, marginals(μ), z)
end

# Marginals with scalar variates transport in a single broadcast:
function transport_from_std(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, z::AbstractVector) where {S<:StdMeasure,M}
    _array_product_from_std(S, μ, z, mspace_flatsize(M))
end
function _array_product_from_std(::Type{S}, μ, z::AbstractVector, ::Tuple{}) where {S}
    mar = marginals(μ)
    broadcast(_FromStd{S}(), mar, maybestatic_reshape(z, maybestatic_size(mar)))
end
function _array_product_from_std(::Type{S}, μ, z::AbstractVector, ::Any) where {S}
    ys, z_rest = _marginals_from_std_with_rest(S, marginals(μ), z)
    if !isempty(z_rest)
        throw(ArgumentError("Length of standard variate doesn't match degrees of freedom of product measure"))
    end
    return ys
end

function _marginals_from_std_with_rest(::Type{S}, νs::Tuple{Vararg{Any}}, z::AbstractVector) where {S}
    y1, z_rest = transport_from_std_with_rest(S, νs[1], z)
    y2_end, z_final_rest = _marginals_from_std_with_rest(S, Base.tail(νs), z_rest)
    return (y1, y2_end...), z_final_rest
end

_marginals_from_std_with_rest(::Type{S}, ::Tuple{}, z::AbstractVector) where {S} = (), z

function _marginals_from_std_with_rest(::Type{S}, νs::AbstractArray{M}, z::AbstractVector) where {S,M}
    idxs = eachindex(νs)
    if isconcretetype(M)
        # The variate type is uniform, so the loop is type stable (the type
        # of the remaining stream stays invariant under repeated view-taking):
        y1, z_rest = transport_from_std_with_rest(S, νs[first(idxs)], z)
        ys = similar(νs, typeof(y1))
        ys[first(idxs)] = y1
        for i in Iterators.drop(idxs, 1)
            ys[i], z_rest = transport_from_std_with_rest(S, νs[i], z_rest)
        end
        return ys, z_rest
    else
        ys_any = Vector{Any}(undef, length(idxs))
        z_rest = z
        for (j, i) in enumerate(idxs)
            ys_any[j], z_rest = transport_from_std_with_rest(S, νs[i], z_rest)
        end
        return [y for y in ys_any], z_rest
    end
end

# Batched transport of array products with scalar-variate marginals in one
# broadcast, the marginals align with the leading dimension of the batch:

function batched_transport_to_std(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, X::AbstractArray) where {S<:StdMeasure,M}
    _array_product_batched_to_std(S, μ, X, mspace_flatsize(M), Val(isconcretetype(M)))
end
function _array_product_batched_to_std(::Type{S}, μ, X::AbstractArray, ::Tuple{}, ::Val{true}) where {S}
    _check_flatsize(X, maybestatic_size(marginals(μ)))
    _as_stream_batch(broadcast(_ToStd{S}(), marginals(μ), X), maybestatic_size(marginals(μ)))
end
function _array_product_batched_to_std(::Type{S}, μ, X::AbstractArray, ::Any, ::Val) where {S}
    _batched_to_std(S, μ, X, mspace_flatsize(μ))
end

function batched_transport_from_std(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, Z::AbstractArray) where {S<:StdMeasure,M}
    _array_product_batched_from_std(S, μ, Z, mspace_flatsize(M), Val(isconcretetype(M)))
end
function _array_product_batched_from_std(::Type{S}, μ, Z::AbstractArray, ::Tuple{}, ::Val{true}) where {S}
    mar = marginals(μ)
    broadcast(_FromStd{S}(), mar, reshape(Z, (map(dynamic, maybestatic_size(mar))..., Base.tail(size(Z))...)))
end
function _array_product_batched_from_std(::Type{S}, μ, Z::AbstractArray, ::Any, ::Val) where {S}
    _batched_from_std(S, μ, Z, mspace_flatsize(μ))
end
