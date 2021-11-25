# Lebesgue measure

export Lebesgue

struct Lebesgue{T} <: PrimitiveMeasure
    support::T
end


sampletype(::Lebesgue) = Float64

testvalue(d::Lebesgue) = testvalue(d.support)

logdensity_def(::Lebesgue, x) = zero(x)

Base.:∘(::typeof(basemeasure), ::Type{Lebesgue}) = Lebesgue

Base.show(io::IO, d::Lebesgue) = print("Lebesgue(",d.support,")")