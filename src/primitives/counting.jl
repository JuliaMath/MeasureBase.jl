# Counting measure

export Counting, CountingBase

struct CountingBase <: PrimitiveMeasure end

@inline mspace_elsize(::CountingBase) = ()
@inline mspace_flatsize(::CountingBase) = ()
@inline mspace_flatsize(::Type{CountingBase}) = ()

insupportof(::CountingBase, x) = true

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
massof(μ::Counting{<:BoundedInts}) = max(μ.support.upper - μ.support.lower + static(1), static(0))
massof(μ::Counting{<:Union{AbstractSet,AbstractRange}}) = length(μ.support)

massof(c::Counting, s::Set) = massof(CountingBase(), filter(insupportof(c), s))

massof(::CountingBase, s::Set) = length(s)

@inline mdomain(μ::Counting{<:SetLike}) = μ.support

# The normalized law on a finite set is the uniform distribution on it,
# drawn by uniform index, the values stay exact:
const _FiniteCounting = Counting{<:Union{AbstractRange,BoundedInts,AbstractSet}}
function batched_rand_impl(ctx::GenContext, μ::_FiniteCounting, sz::SizeLike)
    vals = _counting_values(μ.support)
    n = length(vals)
    n > 0 || throw(ArgumentError("Can't draw random variates from counting measure on an empty set"))
    return _values_at(vals, min.(floor.(Int, _rand_bulk(ctx, sz) .* n), n - 1) .+ 1)
end
_counting_values(r::AbstractRange) = r
_counting_values(s::AbstractSet) = collect(s)
function _counting_values(b::BoundedInts)
    isfinite(b.lower) && isfinite(b.upper) ||
        throw(ArgumentError("Random variates of counting measure require a finite support, got $b"))
    Int(dynamic(b.lower)):Int(dynamic(b.upper))
end
@inline _values_at(vals, i::Integer) = vals[i]
@inline _values_at(vals, idxs::AbstractArray) = map(i -> vals[i], idxs)
