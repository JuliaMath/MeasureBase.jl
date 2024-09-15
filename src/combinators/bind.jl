"""
    struct MeasureBase.Bind{M,K} <: AbstractMeasure

Represents a monatic bind. User code should create instances of `Bind`
directly, but should call `mbind(k, μ)` instead.
"""
struct Bind{M,K} <: AbstractMeasure
    k::K
    μ::M
end

getdof(d::Bind) = NoDOF{typeof(d)}()

function Base.rand(rng::AbstractRNG, ::Type{T}, d::Bind) where {T}
    x = rand(rng, T, d.μ)
    y = rand(rng, T, d.k(x))
    return y
end


"""
    mbind(k, μ)::AbstractMeasure
    
Given

- a measure μ
- a kernel function k that takes values from the support of μ and returns a
  measure

The *monadic bind* operation `mbind(k, μ)` returns is a new measure.

A monadic bind ofen written as `>>=` (e.g. in Haskell), but this symbol is
unavailable in Julia.

```
μ = StdExponential()
ν = mbind(μ) do scale
    pushfwd(Base.Fix1(*, scale), StdNormal())
end
```
"""
mbind(k, μ) = Bind(k, μ)
export mbind
