# Lebesgue measure

export Lebesgue

struct Lebesgue{X} <: PrimitiveMeasure end

Lebesgue(X) = Lebesgue{X}()

gentype(::Lebesgue{ℝ}) = Float64
gentype(::Lebesgue{ℝ₊}) = Float64
gentype(::Lebesgue{𝕀}) = Float64

testvalue(::Lebesgue{ℝ}) = 0.0
testvalue(::Lebesgue{𝕀}) = 0.5
testvalue(::Lebesgue{ℝ₊}) = 1.0
testvalue(::Lebesgue{<:Real}) = 0.0

logdensity_def(::Lebesgue, x) = zero(x)

Base.:∘(::typeof(basemeasure), ::Type{Lebesgue}) = Lebesgue


struct LebesgueCodimOne{D,T,O} <: AbstractMeasure
    ndims ::D
    ortho ::O
end

function logdensityof(d::Density{L1, L2}, x) where {L1<:LebesgueCodimOne, L2<:LebesgueCodimOne}
    μ = d.μ
    ν = d.base
    μ.ndims == ν.ndims || return NaN
    rank([μ.ortho ν.ortho]) == 1 || return NaN
    return 0.0
end

struct LebesgueSimplex <: AbstractMeasure end
