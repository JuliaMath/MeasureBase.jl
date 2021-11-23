# Lebesgue measure

export Lebesgue

struct Lebesgue{X} <: PrimitiveMeasure end

Lebesgue(X) = Lebesgue{X}()

sampletype(::Lebesgue{ℝ}) = Float64
sampletype(::Lebesgue{ℝ₊}) = Float64
sampletype(::Lebesgue{𝕀}) = Float64

testvalue(::Lebesgue{ℝ}) = 0.0
testvalue(::Lebesgue{𝕀}) = 0.5
testvalue(::Lebesgue{ℝ₊}) = 1.0
testvalue(::Lebesgue{<:Real}) = 0.0

logdensity_def(::Lebesgue, x) = zero(x)

Base.:∘(::typeof(basemeasure), ::Type{Lebesgue}) = Lebesgue
