# Lebesgue measure

export Lebesgue

struct LebesgueMeasure <: PrimitiveMeasure end

testvalue(::LebesgueMeasure) = 0.0

struct Lebesgue{T} <: AbstractMeasure
    support::T
end

function Pretty.tile(μ::Lebesgue)
    Pretty.list_layout([Pretty.tile(μ.support)]; prefix=:Lebesgue)
end

gentype(::Lebesgue) = Float64

Lebesgue() = Lebesgue(ℝ)

# basemeasure(::Lebesgue) = LebesgueMeasure()

testvalue(d::Lebesgue) = testvalue(d.support)

proxy(d::Lebesgue) = restrict(in(d.support), LebesgueMeasure())
@useproxy Lebesgue

Base.:∘(::typeof(basemeasure), ::Type{Lebesgue}) = LebesgueMeasure()

Base.show(io::IO, d::Lebesgue) = print(io, "Lebesgue(", d.support, ")")

insupport(μ::Lebesgue, x) = x ∈ μ.support

insupport(::Lebesgue{RealNumbers}, ::Real) = static(true)
