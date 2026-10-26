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

@inline mspace_elsize(μ::Counting) = _valueset_elsize(μ.support)
@inline mspace_flatsize(μ::Counting) = _valueset_flatsize(μ.support)
@inline mspace_flatsize(::Type{<:Counting{IntegerValues}}) = ()
@inline mspace_flatsize(::Type{<:Counting{<:BoundedInts}}) = ()

testvalue(::Type{T}, d::Counting) where {T} = testvalue(T, d.support)

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
