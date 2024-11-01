export AbstractLikelihood, Likelihood

abstract type AbstractLikelihood end

# @inline function logdensityof(ℓ::AbstractLikelihood, p)
#     t() = dynamic(unsafe_logdensityof(ℓ, p))
#     f() = -Inf
#     ifelse(insupport(ℓ, p), t, f)()
# end

# insupport(ℓ::AbstractLikelihood, p) = insupport(_eval_k(ℓ, p), ℓ.x)

@doc raw"""
    Likelihood(k, x)

Default result of [`likelihoodof(k, x)`](@ref).

Likelihoods are most commonly used in conjunction with an existing _prior_
measure to yield a new measure, the _posterior_. In Bayes's Law, we have

``P(θ|x) ∝ P(θ) P(x|θ)``

Here ``P(θ)`` is the prior. If we consider ``P(x|θ)`` as a function on ``θ``,
then it is called a likelihood.

Since measures are most commonly manipulated using `density` and `logdensity`,
it's awkward to commit a (log-)likelihood to using one or the other. To evaluate
a `Likelihood`, we therefore use `density` or `logdensity`, depending on the
circumstances. In the latter case, it is of course acting as a log-density.

For example,

    julia> ℓ = Likelihood(Normal{(:μ,)}, 2.0)
    Likelihood(Normal{(:μ,), T} where T, 2.0)

    julia> density_def(ℓ, (μ=2.0,))
    1.0

    julia> logdensity_def(ℓ, (μ=2.0,))
    -0.0

If, as above, the measure includes the parameter information, we can optionally
leave it out of the second argument in the call to `density` or `logdensity`. 

    julia> density_def(ℓ, 2.0)
    1.0

    julia> logdensity_def(ℓ, 2.0)
    -0.0

With several parameters, things work as expected:
    
    julia> ℓ = Likelihood(Normal{(:μ,:σ)}, 2.0)
    Likelihood(Normal{(:μ, :σ), T} where T, 2.0)
    
    julia> logdensity_def(ℓ, (μ=2, σ=3))
    -1.0986122886681098
    
    julia> logdensity_def(ℓ, (2,3))
    -1.0986122886681098
    
    julia> logdensity_def(ℓ, [2, 3])
    -1.0986122886681098

---------

Finally, let's return to the expression for Bayes's Law, 

``P(θ|x) ∝ P(x|θ) P(θ)``

In measure theory, the product on the right side is the Lebesgue integral
of the likelihood with respect to the prior.

For example, say we have

    μ ~ Normal()
    x ~ Normal(μ,σ)
    σ = 1

and we observe `x=3`. We can compute the (non-normalized) posterior measure on
`μ` as

    julia> prior = Normal()
    julia> likelihood = Likelihood(μ -> Normal(μ, 1), 3)
    julia> post = mintegrate(likelihood, prior)
    julia> post isa MeasureBase.DensityMeasure
    true
    julia> logdensity_rel(post, Lebesgue(), 2)
    -4.337877066409345
"""
struct Likelihood{K,X} <: AbstractLikelihood
    k::K
    x::X

    Likelihood{K,X}(k, x) where {K,X} = new{K,X}(k, x)
end

# For type stability, in case k is a type (resp. a constructor):
Likelihood(k, x::X) where {X} = Likelihood{Core.Typeof(k),X}(k, x)

(lik::AbstractLikelihood)(p) = exp(ULogarithmic, logdensityof(_eval_k(lik, p), lik.x))

_eval_k(ℓ::AbstractLikelihood, p) = asmeasure(ℓ.k(p))

DensityInterface.DensityKind(::AbstractLikelihood) = IsDensity()

function Pretty.quoteof(ℓ::Likelihood)
    k = Pretty.quoteof(ℓ.k)
    x = Pretty.quoteof(ℓ.x)
    :(Likelihood($k, $x))
end

function Base.show(io::IO, ℓ::Likelihood)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, ℓ)
end

insupport(ℓ::AbstractLikelihood, p) = insupport(_eval_k(ℓ, p), ℓ.x)

@inline function logdensityof(ℓ::AbstractLikelihood, p)
    logdensityof(_eval_k(ℓ, p), ℓ.x)
end

@inline function unsafe_logdensityof(ℓ::AbstractLikelihood, p)
    return unsafe_logdensityof(_eval_k(ℓ, p), ℓ.x)
end

# basemeasure(ℓ::Likelihood) = @error "Likelihood requires local base measure"

export likelihoodof

@doc raw"""
    likelihoodof(k, x)

Returns the likelihood of observing `x` under a family of probability
measures that is generated by a transition kernel `k(θ)`.

`k(θ)` maps points in the parameter space to measures (resp. objects that can
be converted to measures) on a implicit set `Χ` that contains values like `x`.

`likelihoodof(k, x)` returns a likelihood object. A likelihhood is **not** a
measure, it is a function from the parameter space to `ℝ₊`. Likelihood
objects can also be interpreted as "generic densities" (but **not** as
probability densities).

`likelihoodof(k, x)` implicitly chooses `ξ  = rootmeasure(k(θ))` as the
reference measure on the observation set `Χ`. Note that this implicit
`ξ` **must** be independent of `θ`.

`ℒₓ = likelihoodof(k, x)` has the mathematical interpretation

```math
\mathcal{L}_x(\theta) = \frac{\rm{d}\, k(\theta)}{\rm{d}\, \chi}(x)
```

`likelihoodof` must return an object that implements the
[`DensityInterface`](https://github.com/JuliaMath/DensityInterface.jl)` API
and `ℒₓ = likelihoodof(k, x)` must satisfy

```julia
log(ℒₓ(θ)) == logdensityof(ℒₓ, θ) ≈ logdensityof(k(θ), x)

DensityKind(ℒₓ) isa IsDensity
```

By default, an instance of [`MeasureBase.Likelihood`](@ref) is returned.
"""
function likelihoodof end

likelihoodof(k, x) = Likelihood(k, x)

###############################################################################
# At the least, we need to think through in some more detail whether
# (log-)likelihood ratios expressed in this way are correct and useful. For now
# this code is commented out; we may remove it entirely in the future.

# export log_likelihood_ratio

# """
#     log_likelihood_ratio(ℓ::Likelihood, p, q)

# Compute the log of the likelihood ratio, in order to compare two choices for
# parameters. This is computed as

#     logdensity_rel(_eval_k(ℓ, p), ℓ.k(q), ℓ.x)

# Since `logdensity_rel` can leave common base measure unevaluated, this can be
# more efficient than

#     logdensityof(_eval_k(ℓ, p), ℓ.x) - logdensityof(ℓ.k(q), ℓ.x)
# """
# log_likelihood_ratio(ℓ::Likelihood, p, q) = logdensity_rel(_eval_k(ℓ, p), ℓ.k(q), ℓ.x)

# # likelihoodof(k, x; kwargs...) = likelihoodof(k, x, NamedTuple(kwargs))

# export likelihood_ratio

# """
#     likelihood_ratio(ℓ::Likelihood, p, q)

# Compute the log of the likelihood ratio, in order to compare two choices for
# parameters. This is equal to

#     density_rel(_eval_k(ℓ, p), ℓ.k(q), ℓ.x)

# but is computed using LogarithmicNumbers.jl to avoid underflow and overflow.
# Since `density_rel` can leave common base measure unevaluated, this can be
# more efficient than

#     logdensityof(_eval_k(ℓ, p), ℓ.x) - logdensityof(ℓ.k(q), ℓ.x)
# """
# function likelihood_ratio(ℓ::Likelihood, p, q)
#     exp(ULogarithmic, logdensity_rel(_eval_k(ℓ, p), ℓ.k(q), ℓ.x))
# end
