# Counting measure

export Counting, CountingMeasure

struct CountingMeasure <: PrimitiveMeasure end

struct Counting{T} <: AbstractMeasure
    support::T
end

Counting() = Counting(ℤ)

basemeasure_type(::Type{C}) where {C<:Counting}= CountingMeasure

testvalue(d::Counting) = testvalue(d.support)

proxy(d::Counting) = restrict(in(d.support), CountingMeasure())

Base.:∘(::typeof(basemeasure), ::Type{Counting}) = CountingMeasure()

Base.show(io::IO, d::Counting) = print(io, "Counting(",d.support,")")

insupport(μ::Counting, x) = x ∈ μ.support

