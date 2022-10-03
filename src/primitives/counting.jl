# Counting measure

export Counting, CountingBase

struct CountingBase <: PrimitiveMeasure end

insupport(::CountingBase, x) = true

struct Counting{T} <: AbstractMeasure
    support::T

    Counting(supp) = new{Core.Typeof(supp)}(supp)
end

function logdensity_def(μ::Counting, x)
    insupport(μ, x) ? 0.0 : -Inf
end

basemeasure(::Counting) = CountingBase()

Counting() = Counting(ℤ)

testvalue(::Type{T}, d::Counting) where {T} = testvalue(T, d.support)

proxy(d::Counting) = restrict(in(d.support), CountingBase())

Base.:∘(::typeof(basemeasure), ::Type{Counting}) = CountingBase()

Base.show(io::IO, d::Counting) = print(io, "Counting(", d.support, ")")

insupport(μ::Counting, x) = x ∈ μ.support

insupport(μ::Counting{T}, x) where {T<:Type} = x isa μ.support

massof(c::Counting, s::Set) = massof(CountingBase(), filter(insupport(c), s))

massof(::CountingBase, s::Set) = length(s)
