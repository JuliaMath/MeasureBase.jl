# Lebesgue measure

export Lebesgue

struct LebesgueMeasure <: PrimitiveMeasure end

testvalue(::LebesgueMeasure) = 0.0


struct Lebesgue{T} <: AbstractMeasure
    support::T
end

Lebesgue() = Lebesgue(ℝ)

gentype(::Lebesgue) = Float64

testvalue(d::Lebesgue) = testvalue(d.support)

proxy(d::Lebesgue) = restrict(in(d.support), LebesgueMeasure())

Base.:∘(::typeof(basemeasure), ::Type{Lebesgue}) = LebesgueMeasure()

Base.show(io::IO, d::Lebesgue) = print("Lebesgue(",d.support,")")

insupport(μ::Lebesgue, x) = x ∈ μ.support

insupport(::Lebesgue{RealNumbers}, ::Real) = static(true)