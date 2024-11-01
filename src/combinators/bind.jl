struct Bind{M,K} <: AbstractMeasure
    μ::M
    k::K
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
mbind(k, μ) = Bind(μ, k)
export mbind

function Base.rand(rng::AbstractRNG, ::Type{T}, d::Bind) where {T}
    x = rand(rng, T, d.μ)
    y = rand(rng, T, d.k(x))
    return y
end


# ToDo: Remove `bind` (breaking).
@noinline function bind(μ, k)
    Base.depwarn("`foo(μ, k)` is deprecated, use `mbind(k, μ)` instead.", :bind)
    mbind(k, μ)
end


# ToDo: Remove `↣` (breaking): It looks too similar to the `>=>` "fish"
# operator (e.g. in Haskell) that is typically understood to take two monadic
# functions as an argument, while a bind take a monad and a monadic functions.
@deprecate ↣(μ, k) mbind(μ, k)
export ↣
