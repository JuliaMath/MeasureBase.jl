# Counting measure

export Counting, CountingBase

"""
    CountingBase()

Counting measure on all of the real numbers, the primitive measure behind
[`Counting`](@ref). It isn't σ-finite, so measures over counting measures
use `Counting(S)` as their base measure.
"""
struct CountingBase <: PrimitiveMeasure end

@inline mspace_elsize(::CountingBase) = ()
@inline mspace_flatsize(::CountingBase) = ()
@inline mspace_flatsize(::Type{CountingBase}) = ()

insupportof(::CountingBase, x) = true

"""
    Counting(S = ℤ)

Counting measure on the set `S`, e.g. a range, or on the values of `S` if
`S` is a type. `Counting(S)` is its own base measure.
"""
struct Counting{T} <: AbstractMeasure
    support::T

    Counting(supp) = new{Core.Typeof(supp)}(supp)
end

@inline function logdensityof_impl(μ::Counting, x::Number)
    R = float(typeof(x))
    _checksupport(insupportof(μ, x), zero(R))
end

@inline logdensityof_impl(μ::Counting, x) = _checksupport(insupportof(μ, x), 0.0)

# Counting measures on a set are their own base measures, so that the base
# carries the set (`CountingBase()` has an atom at every real number and is
# not σ-finite):
basemeasure(μ::Counting) = μ

@inline logdensity_def(::Counting, x) = zero(_logd_numtype(x))

Counting() = Counting(ℤ)

# Supports of scalar variates given as collections or types of numbers:
const _ScalarCountingSupport = Union{AbstractRange{<:Number},AbstractSet{<:Number},AbstractVector{<:Number},Type{<:Number}}
const _ScalarCounting = Counting{<:Union{IntegerValues,BoundedInts,_ScalarCountingSupport}}

@inline mspace_elsize(μ::Counting) = _valueset_elsize(μ.support)
@inline mspace_flatsize(μ::Counting) = _valueset_flatsize(μ.support)
@inline mspace_elsize(::_ScalarCounting) = ()
@inline mspace_flatsize(::_ScalarCounting) = ()
@inline mspace_flatsize(::Type{<:_ScalarCounting}) = ()

testvalue(::Type{T}, d::Counting) where {T} = testvalue(T, d.support)
function testvalue(::Type{T}, d::Counting{<:Union{AbstractRange{<:Number},AbstractSet{<:Number},AbstractVector{<:Number}}}) where {T}
    isempty(d.support) ? zero(T) : convert(T, first(d.support))
end
testvalue(::Type{T}, ::Counting{Type{S}}) where {T,S<:Number} = isconcretetype(S) ? zero(S) : convert(S, zero(T))

proxy(d::Counting) = restrict(in(d.support), CountingBase())

Base.show(io::IO, d::Counting) = print(io, "Counting(", d.support, ")")

insupportof(μ::Counting, x) = x ∈ μ.support

insupportof(μ::Counting{T}, x) where {T<:Type} = x isa μ.support

massof(::Counting{IntegerValues}) = static(Inf)
massof(μ::Counting{<:BoundedInts}) = _integer_cardinality(minimum(μ.support), maximum(μ.support))
function _integer_cardinality(lo, hi)
    count = _integer_count(lo, hi)
    ifelse(_empty_integer_bounds(lo, hi), zero(count), count)
end
_integer_count(lo, hi) = hi - lo + static(1)
_integer_count(lo::Integer, hi::Integer) = widen(hi) - widen(lo) + 1
massof(μ::Counting{<:Union{AbstractSet,AbstractRange}}) = _counting_size(μ.support)
_counting_size(s) = length(s)
_counting_size(s::AbstractUnitRange{<:Integer}) = _integer_cardinality(first(s), last(s))

massof(c::Counting, s::Set) = massof(CountingBase(), filter(insupportof(c), s))

massof(::CountingBase, s::Set) = length(s)

@inline mdomain(μ::Counting{<:SetLike}) = μ.support

# The normalized law on a finite set is the uniform distribution on it,
# drawn by uniform index, the values stay exact:
const _FiniteCounting = Counting{<:Union{AbstractRange,BoundedInts,AbstractSet}}
function batched_rand_impl(ctx::GenContext, μ::_FiniteCounting, sz::SizeLike)
    vals = _counting_values(μ.support)
    isempty(vals) && throw(ArgumentError("Can't draw random variates from counting measure on an empty set"))
    if get_rng(ctx) isa ConstantRNG
        x = _counting_midpoint(vals)
        return isempty(sz) ? x : _const_bulk(ctx, x, sz)
    elseif _on_host(ctx)
        return _host_counting_draw(get_rng(ctx), vals, sz)
    else
        return _device_counting_draw(ctx, vals, sz)
    end
end
_counting_values(r::AbstractRange) = r
_counting_values(s::AbstractSet) = collect(s)
function _counting_values(b::BoundedInts)
    lo, hi = minimum(b), maximum(b)
    isfinite(lo) && isfinite(hi) ||
        throw(ArgumentError("Random variates of counting measure require a finite support, got $b"))
    _counting_endpoint(dynamic(lo)):_counting_endpoint(dynamic(hi))
end
_counting_endpoint(x::Integer) = x
_counting_endpoint(x) = Int(x)
_counting_midpoint(vals) = vals[length(vals) ÷ 2 + 1]
_counting_midpoint(vals::AbstractUnitRange{T}) where {T<:Integer} =
    convert(T, first(vals) + _counting_size(vals) ÷ 2)
_host_counting_draw(rng, vals, ::Tuple{}) = rand(rng, vals)
_host_counting_draw(rng, vals, sz::SizeLike) = rand(rng, vals, asnonstatic(sz))
function _host_counting_draw(rng, vals, sz::StaticSizeLike)
    A = staticarray_type(eltype(vals), canonical_size(sz))
    A(ntuple(_ -> rand(rng, vals), Val(length(A))))
end

function _device_counting_draw(ctx, vals, sz)
    count = _counting_size(vals)
    count <= (UInt128(1) << 64) || throw(ArgumentError("Device counting supports at most 2^64 atoms"))
    n = count % UInt64  # Zero represents the full 64-bit range.
    bits = rand!(ctx, allocate_array(ctx, UInt64, isempty(sz) ? (1,) : asnonstatic(sz)))
    # [threshold, 2^64) contains equally many representatives of each residue.
    if !iszero(n)
        threshold = -n % n
        iszero(threshold) || _reject_counting_bits!(ctx, bits, threshold)
    end
    support = _on_compute_unit(ctx, vals)
    samples = map(x -> _counting_value(support, iszero(n) ? x : x % n), bits)
    return isempty(sz) ? _counting_scalar(samples) : samples
end

_counting_scalar(samples) = only(Array(samples))

_counting_value(vals, offset) = vals[signed(offset) + 1]
# Index arithmetic on full-width ranges wraps in the type of the atoms.
_counting_value(vals::AbstractUnitRange{Int64}, offset) = first(vals) + signed(offset)
_counting_value(vals::AbstractUnitRange{UInt64}, offset) = first(vals) + offset

function _reject_counting_bits!(ctx, bits, threshold)
    while any(<(threshold), bits)
        _redraw_counting_bits!(ctx, bits, threshold)
    end
    return bits
end
function _redraw_counting_bits!(ctx, bits, threshold)
    fresh = rand!(ctx, similar(bits))
    bits .= ifelse.(bits .< threshold, fresh, bits)
end
