export AbstractLikelihood, Likelihood

abstract type AbstractLikelihood end

# @inline function logdensityof(ℓ::AbstractLikelihood, p)
#     t() = dynamic(unsafe_logdensityof(ℓ, p))
#     f() = -Inf
#     ifelse(insupport(ℓ, p), t, f)()
# end

# insupport(ℓ::AbstractLikelihood, p) = insupport(ℓ.k(p), ℓ.x)

@doc raw"""
    Likelihood(k::AbstractTransitionKernel, x)

"Observe" a value `x`, yielding a function from the parameters to ℝ.

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

    Likelihood(M<:ParameterizedMeasure, constraint::NamedTuple, x)

In some cases the measure might have several parameters, and we may want the
(log-)likelihood with respect to some subset of them. In this case, we can use
the three-argument form, where the second argument is a constraint. For example,

    julia> ℓ = Likelihood(Normal{(:μ,:σ)}, (σ=3.0,), 2.0)
    Likelihood(Normal{(:μ, :σ), T} where T, (σ = 3.0,), 2.0)

Similarly to the above, we have

    julia> density_def(ℓ, (μ=2.0,))
    0.3333333333333333

    julia> logdensity_def(ℓ, (μ=2.0,))
    -1.0986122886681098

    julia> density_def(ℓ, 2.0)
    0.3333333333333333

    julia> logdensity_def(ℓ, 2.0)
    -1.0986122886681098

-----------------------

Finally, let's return to the expression for Bayes's Law, 

``P(θ|x) ∝ P(θ) P(x|θ)``

The product on the right side is computed pointwise. To work with this in
MeasureBase, we have a "pointwise product" `⊙`, which takes a measure and a
likelihood, and returns a new measure, that is, the unnormalized posterior that
has density ``P(θ) P(x|θ)`` with respect to the base measure of the prior.

For example, say we have

    μ ~ Normal()
    x ~ Normal(μ,σ)
    σ = 1

and we observe `x=3`. We can compute the posterior measure on `μ` as

    julia> post = Normal() ⊙ Likelihood(Normal{(:μ, :σ)}, (σ=1,), 3)
    Normal() ⊙ Likelihood(Normal{(:μ, :σ), T} where T, (σ = 1,), 3)

    julia> logdensity_def(post, 2)
    -2.5
"""
struct Likelihood{K,X} <: AbstractLikelihood
    k::K
    x::X

    Likelihood(k::K, x::X) where {K<:AbstractTransitionKernel,X} = new{K,X}(k, x)
    Likelihood(k::K, x::X) where {K<:Function,X} = new{K,X}(k, x)
    Likelihood(μ, x) = Likelihood(kernel(μ), x)
end

(lik::AbstractLikelihood)(p) = exp(ULogarithmic, logdensityof(lik.k(p), lik.x))

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

insupport(ℓ::AbstractLikelihood, p) = insupport(ℓ.k(p), ℓ.x)

@inline function logdensityof(ℓ::AbstractLikelihood, p)
    logdensityof(ℓ.k(p), ℓ.x)
end

@inline function unsafe_logdensityof(ℓ::AbstractLikelihood, p)
    return unsafe_logdensityof(ℓ.k(p), ℓ.x)
end

# basemeasure(ℓ::Likelihood) = @error "Likelihood requires local base measure"

export likelihoodof

"""
    likelihoodof(k::AbstractTransitionKernel, x; constraints...)
    likelihoodof(k::AbstractTransitionKernel, x, constraints::NamedTuple)

A likelihood is *not* a measure. Rather, a likelihood acts on a measure, through
the "pointwise product" `⊙`, yielding another measure.
"""
function likelihoodof end

likelihoodof(k, x, ::NamedTuple{()}) = Likelihood(k, x)

likelihoodof(k, x; kwargs...) = likelihoodof(k, x, NamedTuple(kwargs))

likelihoodof(k, x, pars::NamedTuple) = likelihoodof(kernel(k, pars), x)

likelihoodof(k::AbstractTransitionKernel, x) = Likelihood(k, x)

export log_likelihood_ratio

"""
    log_likelihood_ratio(ℓ::Likelihood, p, q)

Compute the log of the likelihood ratio, in order to compare two choices for
parameters. This is computed as

    logdensity_rel(ℓ.k(p), ℓ.k(q), ℓ.x)

Since `logdensity_rel` can leave common base measure unevaluated, this can be
more efficient than

    logdensityof(ℓ.k(p), ℓ.x) - logdensityof(ℓ.k(q), ℓ.x)
"""
log_likelihood_ratio(ℓ::Likelihood, p, q) = logdensity_rel(ℓ.k(p), ℓ.k(q), ℓ.x)

# likelihoodof(k, x; kwargs...) = likelihoodof(k, x, NamedTuple(kwargs))

export likelihood_ratio

"""
    likelihood_ratio(ℓ::Likelihood, p, q)

Compute the log of the likelihood ratio, in order to compare two choices for
parameters. This is equal to

    density_rel(ℓ.k(p), ℓ.k(q), ℓ.x)

but is computed using LogarithmicNumbers.jl to avoid underflow and overflow.
Since `density_rel` can leave common base measure unevaluated, this can be
more efficient than

    logdensityof(ℓ.k(p), ℓ.x) - logdensityof(ℓ.k(q), ℓ.x)
"""
function likelihood_ratio(ℓ::Likelihood, p, q)
    exp(ULogarithmic, logdensity_rel(ℓ.k(p), ℓ.k(q), ℓ.x))
end
