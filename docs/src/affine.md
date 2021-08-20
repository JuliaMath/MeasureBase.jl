# Affine Transformations

It's very common for measures to be parameterized by `μ` and `σ`, for example as in `Normal(μ=3, σ=4)` or `StudentT(ν=1, μ=3, σ=4)`. In this context, `μ` and `σ` do not always refer to the mean and standard deviation (the `StudentT` above is equivalent to a Cauchy, so both are undefined).

Rather, `μ` is a "location parameter", and `σ` is a "scale parameter". Together these determine an affine transformation

```math
f(z) = σ z + μ
```

Here are below, we'll use ``z`` to represent an "un-transformed" variable, typically coming from a measure like `Normal()` with no location or scale parameters.

Affine transforms are often incorrectly referred to as "linear". Linearity requires ``f(ax + by) = a f(x) + b f(y)`` for scalars ``a`` and ``b``, which only holds for the above ``f`` if ``μ=0``.


## Cholesky-based parameterizations

If the "un-transformed" `z` is a scalar, things are relatively simple. But it's important our approach handle the multivariate case as well.

In the literature, it's common for a multivariate normal distribution to be parameterized by a mean `μ` and covariance matrix `Σ`. This is mathematically convenient, but can be very awkward from a computational perspective.

While MeasureTheory.jl includes (or will include) a parameterization using `Σ`, we prefer to work in terms of its Cholesky decomposition ``σ``.

Using "``σ``" for this may seem strange at first, so we should explain the notation. Let ``σ`` be a lower-triangular matrix satisfying

```math
σ σᵗ = Σ
```

Then given a (multivariate) standard normal ``z``, the covariance matrix of ``σ z + μ`` is

```math
𝕍[σ z + μ] = Σ
```

This is similar to the one dimensional case where

```math
𝕍[σ z + μ] = σ² ,
```

and so the lower Cholesky factor of the covariance generalizes the concept of standard deviation, justifying the notation.

## `Affine` and `AffineTransform`



unif = ∫(x -> 0<x<1, Lebesgue(ℝ))
    f = AffineTransform((μ=3,σ=2))
    g = AffineTransform((μ=3,ω=2))

So for example, the implementation of `StudentT(ν=1, μ=3, σ=4)` is equivalent to

```julia
StudentT(nt::NamedTuple{(:ν,:μ,:σ)}) = Affine((μ=nt.μ, σ=nt.σ), StudentT((ν=1)))
```

