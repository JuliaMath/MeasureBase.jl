# Counting measure

export Counting, CountingMeasure

struct CountingMeasure <: PrimitiveMeasure end

struct Counting{T} <: AbstractMeasure
    support::T
end

function logdensity_def(μ::Counting, x)
    insupport(μ, x) ? 0.0 : -Inf
end

basemeasure(::Counting) = CountingMeasure()

Counting() = Counting(ℤ)

testvalue(d::Counting) = testvalue(d.support)

proxy(d::Counting) = restrict(in(d.support), CountingMeasure())

Base.:∘(::typeof(basemeasure), ::Type{Counting}) = CountingMeasure()

Base.show(io::IO, d::Counting) = print(io, "Counting(", d.support, ")")

insupport(μ::Counting, x) = x ∈ μ.support

insupport(μ::Counting{T}, x) where {T<:Type} = x isa T