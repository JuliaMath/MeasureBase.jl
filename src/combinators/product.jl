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
@inline _marginal_rand(ctx::GenContext, d) = convert_realtype(get_precision(ctx), rand(get_rng(ctx), d))

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

# Batches of tuple and named tuple variates are tuples resp. named tuples
# of batches, the marginal densities add up lazily:
for (bhead, head) in [(:batched_logdensityof_impl, :logdensityof_impl), (:batched_logdensity_def, :logdensity_def)]
    @eval @inline function $bhead(μ::ProductMeasure{<:Tuple}, X::Tuple)
        _lazy_sum(map((m, Xi) -> _batched_kernel($head, m, Xi), marginals(μ), X))
    end
    @eval @inline function $bhead(μ::ProductMeasure{<:NamedTuple{names}}, X::NamedTuple{names}) where {names}
        _lazy_sum(map((m, Xi) -> _batched_kernel($head, m, Xi), values(marginals(μ)), values(X)))
    end
end
@inline _lazy_sum(ℓs::Tuple) = reduce(_lazy_add, ℓs)

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

@inline function mspace_ndims(::Type{<:ProductMeasure{<:AbstractArray{M,N}}}) where {M,N}
    _add_ndims(mspace_ndims(M), N)
end
@inline fixed_stream_size(::Type{<:ProductMeasure{<:AbstractArray{M}}}) where {M} = fixed_stream_size(M)

# Batched kernels over flat storage `(marginal variate dims..., product
# dims..., batch dims...)`. Marginals with scalar variates align with the
# leading dimensions of the batch, so one broadcast evaluates all marginal
# densities. Marginals with array variates are evaluated one by one over
# their slices of the batch.
for (bhead, head) in [(:batched_logdensityof_impl, :logdensityof_impl), (:batched_logdensity_def, :logdensity_def)]
    @eval @inline function $bhead(μ::ProductMeasure{<:AbstractArray{M}}, X::AbstractArray) where {M}
        _array_product_kernel($head, μ, X, _static_ndims_of(mspace_ndims(M)))
    end
end

@inline function _array_product_kernel(f::F, μ::ProductMeasure, X::AbstractArray, k::Integer) where {F}
    _array_product_kernel(f, μ, X, static(k))
end
@inline function _array_product_kernel(f::F, μ::ProductMeasure, X::AbstractArray, ::StaticInteger{0}) where {F}
    mar = marginals(μ)
    _check_flatsize(X, maybestatic_size(mar))
    _sum_leading_dims(_marginal_broadcast(_DynamicPointLogd(f), mar, X), static(ndims(mar)))
end
@inline function _array_product_kernel(f::F, μ::ProductMeasure, X::AbstractArray, ::StaticInteger{K}) where {F,K}
    _marginal_slices_ld(f, marginals(μ), X, Val(K), Val(ndims(X) - K - ndims(marginals(μ))))
end
@noinline function _array_product_kernel(::F, μ::ProductMeasure{<:AbstractArray{M}}, ::AbstractArray, ::NoMSpaceElementSize) where {F,M}
    throw(ArgumentError("Batched density evaluation of products over arrays of marginals of type $(nameof(M)) requires MeasureBase.mspace_ndims to be declared for that type"))
end

struct _DynamicPointLogd{F} <: Function
    f::F
end
@inline (k::_DynamicPointLogd)(m, x) = _dynamic_logd(k.f(m, x), x)

function _marginal_slices_ld(f::F, mar::AbstractArray{<:Any,N}, A::AbstractArray, ::Val{K}, ::Val{B}) where {F,N,K,B}
    if ndims(A) != K + N + B || ntuple(i -> size(A, K + i), Val(N)) != size(mar)
        _throw_size_mismatch()
    end
    lead = ntuple(_ -> Colon(), Val(K))
    trail = ntuple(_ -> Colon(), Val(B))
    ld(i) = _materialize(_batched_kernel(f, mar[i], view(A, lead..., Tuple(i)..., trail...)))
    init = _zero_logd(A, ntuple(i -> size(A, K + N + i), Val(B)))
    return mapreduce(ld, +, CartesianIndices(mar); init = init)
end

@inline _zero_logd(A::AbstractArray, ::Tuple{}) = zero(_logd_numtype(A))
@inline _zero_logd(A::AbstractArray, dims::Tuple) = fill!(similar(A, _logd_numtype(A), dims), 0)

# Point densities of array products at numeric variates go through the
# batched kernel where the variate rank of the marginals is known,
# marginal by marginal otherwise:
@inline _point_ld(f::F, μ::AbstractProductMeasure, x::AbstractArray{<:Number}) where {F} = f(μ, x)
@inline function logdensityof_impl(μ::ProductMeasure{<:AbstractArray{M}}, x::AbstractArray{<:Number}) where {M}
    _array_product_ld(logdensityof_impl, μ, x, mspace_ndims(M))
end
@inline function logdensity_def(μ::ProductMeasure{<:AbstractArray{M}}, x::AbstractArray{<:Number}) where {M}
    _array_product_ld(logdensity_def, μ, x, mspace_ndims(M))
end
@inline function _array_product_ld(f::F, μ::ProductMeasure, x::AbstractArray, ::Integer) where {F}
    _point_result(_materialize(_batched_kernel(f, μ, x)), μ)
end
@inline function _array_product_ld(f::F, μ::ProductMeasure, x::AbstractArray, ::NoMSpaceElementSize) where {F}
    _array_product_ld_nested(f, μ, x)
end
@inline function _array_product_ld_nested(f::F, μ::ProductMeasure, x::AbstractArray) where {F}
    _check_marginal_count(marginals(μ), x)
    mapreduce(_PointLogd(f, nothing), +, marginals(μ), x)
end
@inline (k::_PointLogd{F,Nothing})(m, x) where {F} = _point_ld(k.f, m, x)

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
@inline _sum_dofs(f, mar::AbstractArray{M}) where {M} = _sum_dofs(f, mar, _unit_dof(M))
@inline _sum_dofs(f, mar::AbstractArray, ::True) = length(mar)
@inline _sum_dofs(f, mar::AbstractArray, ::False) = mapreduce(_dynamic_dof ∘ f, +, mar; init = 0)

# Marginals with scalar variates and a standard transport have one degree
# of freedom each, so their total needs no reduction over the marginals
# (which may live on a device):
@inline function _unit_dof(::Type{M}) where {M}
    static(mspace_ndims(M) === 0 && preferred_stdmeasure(M) isa Type{<:StdMeasure})
end
@inline _sum_dofs(f, mar::StaticArray) = mapreduce(f, +, mar; init = static(0))
@inline _dynamic_dof(n::IntegerLike) = dynamic(n)
@inline _dynamic_dof(nodof::AbstractNoDOF) = nodof

function checked_arg(μ::ProductMeasure{<:NTuple{N,Any}}, x::NTuple{N,Any}) where {N}
    map(checked_arg, marginals(μ), x)
end

# Variates of array products are arrays of marginal variates or, for
# marginals with array variates of declared rank, their flat storage:
@propagate_inbounds function checked_arg(μ::ProductMeasure{<:AbstractArray{M}}, x::AbstractArray) where {M}
    @boundscheck _check_product_arg(marginals(μ), x, _static_ndims_of(mspace_ndims(M)))
    return x
end

@inline function _check_product_arg(mar, x::AbstractArray, ::Any)
    _check_marginal_count(mar, x)
    foreach(checked_arg, mar, x)
    return nothing
end
@inline function _check_product_arg(mar::AbstractArray{<:Any,N}, x::AbstractArray{<:Number}, ::StaticInteger{0}) where {N}
    _check_marginal_count(mar, x)
end
@inline function _check_product_arg(mar::AbstractArray{<:Any,N}, x::AbstractArray{<:Number}, ::StaticInteger{K}) where {N,K}
    if ndims(x) != K + N || ntuple(i -> size(x, K + i), Val(N)) != size(mar)
        _throw_marginal_mismatch()
    end
    return nothing
end
@inline _check_product_arg(mar, x::AbstractArray{<:Number}, ::NoMSpaceElementSize) = _check_marginal_count(mar, x)


function checked_arg(
    μ::ProductMeasure{<:NamedTuple{names}},
    x::NamedTuple{names},
) where {names}
    NamedTuple{names}(map(checked_arg, values(marginals(μ)), values(x)))
end


# Transport marginal by marginal, the standard variates of the marginals
# are concatenated in order. Batches of tuple and named tuple variates are
# tuples resp. named tuples of batches.

function transport_to_std(::Type{S}, μ::ProductMeasure{<:Tuple}, x::Tuple) where {S<:StdMeasure}
    _flatten_to_rv(map((m, xi) -> _as_stdstream(transport_to_std(S, m, xi)), marginals(μ), x))
end

function transport_to_std(::Type{S}, μ::ProductMeasure{<:NamedTuple{names}}, x::NamedTuple{names}) where {S<:StdMeasure,names}
    transport_to_std(S, productmeasure(values(marginals(μ))), values(x))
end

function batched_transport_to_std(::Type{S}, μ::ProductMeasure{<:Tuple}, X::Tuple) where {S<:StdMeasure}
    _vcat_std(map((m, Xi) -> batched_transport_to_std(S, m, Xi), marginals(μ), X))
end

function batched_transport_to_std(::Type{S}, μ::ProductMeasure{<:NamedTuple{names}}, X::NamedTuple{names}) where {S<:StdMeasure,names}
    batched_transport_to_std(S, productmeasure(values(marginals(μ))), values(X))
end

@inline _vcat_std(Zs::Tuple) = vcat(Zs...)
@inline _vcat_std(::Tuple{}) = SVector{0,Bool}()

function transport_from_std(::Type{S}, μ::ProductMeasure{<:Union{Tuple,NamedTuple}}, z::AbstractVector) where {S<:StdMeasure}
    x, z_rest = transport_from_std_with_rest(S, μ, z)
    isempty(z_rest) || _throw_std_length_mismatch()
    return x
end

function transport_from_std_with_rest(::Type{S}, μ::ProductMeasure{<:Tuple}, z::AbstractVector) where {S<:StdMeasure}
    _marginals_from_std_with_rest(S, marginals(μ), z)
end

function transport_from_std_with_rest(::Type{S}, μ::ProductMeasure{<:NamedTuple{names}}, z::AbstractVector) where {S<:StdMeasure,names}
    ys, z_rest = _marginals_from_std_with_rest(S, values(marginals(μ)), z)
    return NamedTuple{names}(ys), z_rest
end

function batched_transport_from_std(::Type{S}, μ::ProductMeasure{<:Union{Tuple,NamedTuple}}, Z::AbstractArray) where {S<:StdMeasure}
    X, Z_rest = batched_transport_from_std_with_rest(S, μ, Z, ())
    size(Z_rest, 1) == 0 || _throw_std_length_mismatch()
    return X
end

function batched_transport_from_std_with_rest(::Type{S}, μ::ProductMeasure{<:Tuple}, Z::AbstractArray, sz::Dims) where {S<:StdMeasure}
    _tuple_product_from_std_with_rest(S, μ, Z, sz)
end

function batched_transport_from_std_with_rest(::Type{S}, μ::ProductMeasure{<:NamedTuple{names}}, Z::AbstractArray, sz::Dims) where {S<:StdMeasure,names}
    Xs, Z_rest = _tuple_product_from_std_with_rest(S, productmeasure(values(marginals(μ))), Z, sz)
    return NamedTuple{names}(Xs), Z_rest
end

# One variate per stream is consumed marginal by marginal, several per
# stream via the degrees of freedom of the whole product:
function _tuple_product_from_std_with_rest(::Type{S}, μ, Z::AbstractArray, ::Tuple{}) where {S}
    _marginals_batched_from_std_with_rest(S, marginals(μ), Z)
end
function _tuple_product_from_std_with_rest(::Type{S}, μ, Z::AbstractArray, sz::Dims) where {S}
    _batched_from_std_bydof(S, μ, Z, sz, fast_dof(μ))
end

function _marginals_batched_from_std_with_rest(::Type{S}, νs::Tuple{Vararg{Any}}, Z::AbstractArray) where {S}
    X1, Z_rest = batched_transport_from_std_with_rest(S, νs[1], Z, ())
    X2_end, Z_final_rest = _marginals_batched_from_std_with_rest(S, Base.tail(νs), Z_rest)
    return (X1, X2_end...), Z_final_rest
end

_marginals_batched_from_std_with_rest(::Type{S}, ::Tuple{}, Z::AbstractArray) where {S} = (), Z

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
        ys = similar(Array{typeof(y1)}, axes(νs))
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


# Array products: marginals of scalar variates with one degree of freedom
# each transport in a single broadcast (the marginals align with the leading
# dimensions of the batch), other marginals one by one over their slices of
# the batch.

@inline _fused_marginals(::Type{M}) where {M} = Val(isconcretetype(M) && _unit_dof(M) === static(true))

function batched_transport_to_std(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, X::AbstractArray) where {S<:StdMeasure,M}
    _array_product_batched_to_std(S, μ, X, _fused_marginals(M), _static_ndims_of(mspace_ndims(M)))
end
function _array_product_batched_to_std(::Type{S}, μ, X::AbstractArray, ::Val{true}, ::Any) where {S}
    mar = marginals(μ)
    _check_flatsize(X, maybestatic_size(mar))
    _as_stream_batch(_materialize(_marginal_broadcast(_ToStd{S}(), mar, X)), static(ndims(mar)))
end
function _array_product_batched_to_std(::Type{S}, μ, X::AbstractArray, ::Val{false}, ::StaticInteger{K}) where {S,K}
    mar = marginals(μ)
    n_batch = ndims(X) - K - ndims(mar)
    n_batch >= 0 || _throw_size_mismatch()
    _marginals_to_std_loop(S, mar, X, Val(K), Val(n_batch))
end
@noinline function _array_product_batched_to_std(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, ::AbstractArray, ::Val{false}, ::NoMSpaceElementSize) where {S,M}
    throw(ArgumentError("Batched transport of products over arrays of marginals of type $(nameof(M)) requires MeasureBase.mspace_ndims to be declared for that type"))
end

function _marginals_to_std_loop(::Type{S}, mar::AbstractArray{<:Any,N}, X::AbstractArray, ::Val{K}, ::Val{B}) where {S,N,K,B}
    ntuple(i -> size(X, K + i), Val(N)) == size(mar) || _throw_size_mismatch()
    lead = ntuple(_ -> Colon(), Val(K))
    trail = ntuple(_ -> Colon(), Val(B))
    zs = map(i -> batched_transport_to_std(S, mar[i], view(X, lead..., Tuple(i)..., trail...)), vec(CartesianIndices(mar)))
    isempty(zs) ? similar(X, (0, ntuple(i -> size(X, K + N + i), Val(B))...)) : reduce(vcat, zs)
end

function batched_transport_from_std(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, Z::AbstractArray) where {S<:StdMeasure,M}
    _array_product_batched_from_std(S, μ, Z, _fused_marginals(M), _static_ndims_of(mspace_ndims(M)))
end
function _array_product_batched_from_std(::Type{S}, μ, Z::AbstractArray, ::Val{true}, ::Any) where {S}
    mar = marginals(μ)
    size(Z, 1) == length(mar) || _throw_std_length_mismatch()
    _materialize(_marginal_broadcast(_FromStd{S}(), mar, _reshape_batch(Z, (_batch_dims(mar)..., Base.tail(_batch_dims(Z))...))))
end
function _array_product_batched_from_std(::Type{S}, μ, Z::AbstractArray, ::Val{false}, ::StaticInteger{K}) where {S,K}
    X, Z_rest = _marginals_from_std_loop(S, marginals(μ), Z, Val(K))
    size(Z_rest, 1) == 0 || _throw_std_length_mismatch()
    return X
end
@noinline function _array_product_batched_from_std(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, ::AbstractArray, ::Val{false}, ::NoMSpaceElementSize) where {S,M}
    throw(ArgumentError("Batched transport to products over arrays of marginals of type $(nameof(M)) requires MeasureBase.mspace_ndims to be declared for that type"))
end

function batched_transport_from_std_with_rest(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, Z::AbstractArray, sz::Dims) where {S<:StdMeasure,M}
    _array_product_batched_from_std_with_rest(S, μ, Z, sz, _fused_marginals(M), _static_ndims_of(mspace_ndims(M)))
end
function _array_product_batched_from_std_with_rest(::Type{S}, μ, Z::AbstractArray, sz::Dims, ::Val{true}, ::Any) where {S}
    _batched_from_std_bydof(S, μ, Z, sz, length(marginals(μ)))
end
function _array_product_batched_from_std_with_rest(::Type{S}, μ, Z::AbstractArray, ::Tuple{}, ::Val{false}, ::StaticInteger{K}) where {S,K}
    _marginals_from_std_loop(S, marginals(μ), Z, Val(K))
end
function _array_product_batched_from_std_with_rest(::Type{S}, μ, Z::AbstractArray, sz::Dims, ::Val{false}, ::StaticInteger{K}) where {S,K}
    _batched_from_std_bydof(S, μ, Z, sz, fast_dof(μ))
end
@noinline function _array_product_batched_from_std_with_rest(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, ::AbstractArray, ::Dims, ::Val{false}, ::NoMSpaceElementSize) where {S,M}
    throw(ArgumentError("Batched transport to products over arrays of marginals of type $(nameof(M)) requires MeasureBase.mspace_ndims to be declared for that type"))
end

# The marginals consume the streams one after the other, their variates
# fill the batch `(marginal variate dims..., product dims..., batch dims...)`:
function _marginals_from_std_loop(::Type{S}, mar::AbstractArray{<:Any,N}, Z::AbstractArray, ::Val{K}) where {S,N,K}
    idxs = vec(CartesianIndices(mar))
    batch_dims = Base.tail(size(Z))
    lead = ntuple(_ -> Colon(), Val(K))
    trail = ntuple(_ -> Colon(), Val(length(batch_dims)))
    if isempty(idxs)
        return similar(Z, (ntuple(_ -> 0, Val(K))..., size(mar)..., batch_dims...)), Z
    end
    X1, Z_rest = batched_transport_from_std_with_rest(S, mar[idxs[1]], Z, ())
    X = similar(Z, eltype(X1), (ntuple(i -> size(X1, i), Val(K))..., size(mar)..., batch_dims...))
    X[lead..., Tuple(idxs[1])..., trail...] = X1
    for i in idxs[2:end]
        Xi, Z_rest = batched_transport_from_std_with_rest(S, mar[i], Z_rest, ())
        X[lead..., Tuple(i)..., trail...] = Xi
    end
    return X, Z_rest
end

# Point transport of array products: arrays of marginal variates with flat
# storage and flat variates go through the batched kernels, marginals
# without a declared variate rank transport one by one.
function transport_to_std(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, x::AbstractArray) where {S<:StdMeasure,M}
    _array_product_to_std(S, μ, x, _flat_storage(x), _static_ndims_of(mspace_ndims(M)))
end
@inline function _array_product_to_std(::Type{S}, μ, x::AbstractArray, x_flat::AbstractArray, ::StaticInteger) where {S}
    _single_std(batched_transport_to_std(S, μ, x_flat))
end
function _array_product_to_std(::Type{S}, μ, x::AbstractArray, ::Any, ::Any) where {S}
    _check_marginal_count(marginals(μ), x)
    zs = [_as_stdstream(transport_to_std(S, m, xi)) for (m, xi) in zip(marginals(μ), x)]
    isempty(zs) ? SVector{0,Bool}() : reduce(vcat, zs)
end

# Marginals with variates of fixed size and declared rank yield a nested
# view of the flat variate batch, others transport marginal by marginal:
function transport_from_std(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, z::AbstractVector) where {S<:StdMeasure,M}
    _array_product_from_std(S, μ, z, fixed_stream_size(M), _static_ndims_of(mspace_ndims(M)))
end
@inline function _array_product_from_std(::Type{S}, μ, z::AbstractVector, ::True, k::StaticInteger) where {S}
    _nest_leaf(batched_transport_from_std(S, μ, z), k)
end
function _array_product_from_std(::Type{S}, μ, z::AbstractVector, ::Any, ::Any) where {S}
    ys, z_rest = _marginals_from_std_with_rest(S, marginals(μ), z)
    isempty(z_rest) || _throw_std_length_mismatch()
    return ys
end

function transport_from_std_with_rest(::Type{S}, μ::ProductMeasure{<:AbstractArray{M}}, z::AbstractVector) where {S<:StdMeasure,M}
    _array_product_from_std_with_rest(S, μ, z, fixed_stream_size(M), _static_ndims_of(mspace_ndims(M)))
end
function _array_product_from_std_with_rest(::Type{S}, μ, z::AbstractVector, ::True, k::StaticInteger) where {S}
    X, z_rest = batched_transport_from_std_with_rest(S, μ, z, ())
    return _nest_leaf(X, k), z_rest
end
function _array_product_from_std_with_rest(::Type{S}, μ, z::AbstractVector, ::Any, ::Any) where {S}
    _marginals_from_std_with_rest(S, marginals(μ), z)
end


# Streams: tuple products consume marginal by marginal, so marginals of
# value-dependent size are supported for a single variate per stream.
function batched_logdensityof_with_rest(μ::ProductMeasure{<:Tuple}, X::AbstractArray, ::Tuple{})
    _marginals_ld_with_rest(marginals(μ), X)
end
function batched_logdensityof_with_rest(μ::ProductMeasure{<:Tuple}, x::AbstractVector, ::Tuple{})
    _marginals_ld_with_rest(marginals(μ), x)
end
function _marginals_ld_with_rest(ms::Tuple, X::AbstractArray)
    ℓ1, X2 = batched_logdensityof_with_rest(ms[1], X, ())
    ℓ_rest, X_rest = _marginals_ld_with_rest(Base.tail(ms), X2)
    return _lazy_add(ℓ1, ℓ_rest), X_rest
end
function _marginals_ld_with_rest(ms::Tuple{Any}, X::AbstractArray)
    batched_logdensityof_with_rest(ms[1], X, ())
end
function batched_logdensityof_with_rest(μ::ProductMeasure{<:NamedTuple{names}}, X::AbstractArray, sz::Tuple{}) where {names}
    batched_logdensityof_with_rest(productmeasure(values(marginals(μ))), X, sz)
end
function batched_logdensityof_with_rest(μ::ProductMeasure{<:NamedTuple{names}}, x::AbstractVector, sz::Tuple{}) where {names}
    batched_logdensityof_with_rest(productmeasure(values(marginals(μ))), x, sz)
end

@inline function fixed_stream_size(::Type{<:ProductMeasure{M}}) where {M<:Tuple}
    static(all(T -> fixed_stream_size(T) === static(true), M.parameters))
end
@inline function fixed_stream_size(::Type{<:ProductMeasure{NamedTuple{names,M}}}) where {names,M<:Tuple}
    fixed_stream_size(ProductMeasure{M})
end


# Broadcasts over struct arrays of marginals run over their leaf columns,
# the marginals are rebuilt from the column values inside the kernel:

@inline _leaf_columns(sa::StructArray) = _leaf_columns_of(values(StructArrays.components(sa)))
@inline _leaf_columns_of(cs::Tuple) = (_leaf_columns_of(first(cs))..., _leaf_columns_of(Base.tail(cs))...)
@inline _leaf_columns_of(::Tuple{}) = ()
@inline _leaf_columns_of(c::StructArray) = _leaf_columns(c)
@inline _leaf_columns_of(c::AbstractArray{T}) where {T} = Base.issingletontype(T) ? () : (c,)

@generated function _rebuild_element(::Type{SA}, vals::Tuple) where {SA<:StructArray}
    expr, _ = _rebuild_expr(SA, 1)
    return expr
end
function _rebuild_expr(::Type{SA}, i::Int) where {T,N,C,SA<:StructArray{T,N,C}}
    args = Any[]
    for CT in C.parameters[2].parameters
        if CT <: StructArray
            e, i = _rebuild_expr(CT, i)
            push!(args, e)
        elseif Base.issingletontype(eltype(CT))
            push!(args, :($(eltype(CT).instance)))
        else
            push!(args, :(vals[$i]))
            i += 1
        end
    end
    return :(constructorof($T)($(args...))), i
end

struct _WithElement{SA,G} <: Function
    g::G
end
_WithElement{SA}(g::G) where {SA,G} = _WithElement{SA,G}(g)
@inline function (k::_WithElement{SA})(args::Vararg{Any,N}) where {SA,N}
    k.g(_rebuild_element(SA, Base.front(args)), args[end])
end

@inline function _marginal_broadcast(g::G, mar::AbstractArray, X) where {G}
    Broadcast.instantiate(Broadcast.broadcasted(g, mar, X))
end
@inline function _marginal_broadcast(g::G, mar::StructArray, X) where {G}
    Broadcast.instantiate(Broadcast.broadcasted(_WithElement{typeof(mar)}(g), _leaf_columns(mar)..., X))
end

Adapt.adapt_structure(to, μ::ProductMeasure) = ProductMeasure(Adapt.adapt(to, marginals(μ)))
