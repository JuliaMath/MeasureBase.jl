"""
    module MeasureOperators

Defines the following operators for measures:

* `f ⋄ μ == pushfwd(f, μ)`
* `μ ⊙ f == pullbck(f, μ)`
* `μ ▷ k == mbind(k, μ)`
* `⊗(μs...) == productmeasure(μs)`
* `∫(f, μ) == mintegrate(f, μ)`
* `∫exp(f, μ) == mintegrate_exp(f, μ)`
* `𝒹(ν, μ) == density_rel(ν, μ)`
* `log𝒹(ν, μ) == logdensity_rel(ν, μ)`
"""
module MeasureOperators

using MeasureBase: AbstractMeasure
using MeasureBase: pushfwd, pullbck, mbind, productmeasure
using MeasureBase: mintegrate, mintegrate_exp, density_rel, logdensity_rel
using InverseFunctions: inverse

@doc raw"""
    ⋄(f, μ::AbstractMeasure) = pushfwd(f, μ)

The `\\diamond` operator denotes a pushforward operation: `ν = f ⋄ μ`
generates a
[pushforward measure](https://en.wikipedia.org/wiki/Pushforward_measure).

A common mathematical notation for a pushforward is ``f_*μ``, but as
there is no "subscript-star" operator in Julia, we use `⋄`.

See [`pushfwd(f, μ)`](@ref) for details.

Also see [`ν ⊙ f`](@ref), the pullback operator.
"""
⋄(f, μ::AbstractMeasure) = pushfwd(f, μ)
export ⋄

@doc raw"""
    ⊙(ν::AbstractMeasure, f) = pullbck(f, ν)

The `\\odot` operator denotes a pullback operation.

See also [`pullbck(ν, f)`](@ref) for details. Note that `pullbck` takes it's
arguments in different order, in keeping with the Julia convention of
passing functions as the first argument. A pullback is mathematically the
precomposition of a measure `μ`` with the function `f` applied to sets. so
`⊙` takes the measure as the first and the function as the second argument,
as common in mathematical notation for precomposition.

A common mathematical notation for pullback in measure theory is
``f \circ μ``, but as `∘` is used for function composition in Julia and as
`f` semantically acts point-wise on sets, we use `⊙`.

Also see [f ⋄ μ](@ref), the pushforward operator.
"""
⊙(ν::AbstractMeasure, f) = pullbck(f, ν)
export ⊙

"""
    μ ▷ k = mbind(k, μ)

The `\\triangleright` operator denotes a measure monadic bind operation.

A common operator choice for a monadic bind operator is `>>=` (e.g. in
the Haskell programming language), but this has a different meaning in
Julia and there is no close equivalent, so we use `▷`.

See [`mbind(k, μ)`](@ref) for details. Note that `mbind` takes its
arguments in different order, in keeping with the Julia convention of
passing functions as the first argument. `▷`, on the other hand, takes
its arguments in the order common for monadic binds in functional
programming (like the Haskell `>>=` operator) and mathematics.
"""
▷(μ::AbstractMeasure, k) = mbind(k, μ)
export ▷

# ToDo: Use `⨂` instead of `⊗` for better readability?
"""
    ⊗(μs::AbstractMeasure...) = productmeasure(μs)

`⊗` is an operator for building product measures.

See [`productmeasure(μs)`](@ref) for details.
"""
⊗(μs::AbstractMeasure...) = productmeasure(μs)
export ⊗

"""
    ∫(f, μ::AbstractMeasure) = mintegrate(f, μ)

Denotes an indefinite integral of the function `f` with respect to the
measure `μ`.

See [`mintegrate(f, μ)`](@ref) for details.
"""
∫(f, μ::AbstractMeasure) = mintegrate(f, μ)
export ∫

"""
    ∫exp(f, μ::AbstractMeasure) = mintegrate_exp(f, μ)

Generates a new measure that is the indefinite integral of `exp` of `f`
with respect to the measure `μ`.

See [`mintegrate_exp(f, μ)`](@ref) for details.
"""
∫exp(f, μ::AbstractMeasure) = mintegrate_exp(f, μ)
export ∫exp

"""
    𝒹(ν, μ) = density_rel(ν, μ)

Compute the density, i.e. the
[Radom-Nikodym derivative](https://en.wikipedia.org/wiki/Radon%E2%80%93Nikodym_theorem)
of `ν`` with respect to `μ`.

For details, see [`density_rel(ν, μ)`}(@ref).
"""
𝒹(ν, μ::AbstractMeasure) = density_rel(ν, μ)
export 𝒹

"""
    log𝒹(ν, μ) = logdensity_rel(ν, μ)

Compute the log-density, i.e. the logarithm of the
[Radom-Nikodym derivative](https://en.wikipedia.org/wiki/Radon%E2%80%93Nikodym_theorem)
of `ν`` with respect to `μ`.

For details, see [`logdensity_rel(ν, μ)`}(@ref).
"""
log𝒹(ν, μ::AbstractMeasure) = logdensity_rel(ν, μ)
export log𝒹

end # module MeasureOperators
