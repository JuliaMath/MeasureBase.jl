"""
    abstract type AbstractLikelihood <: Function

Abstract supertype for likelihood objects.

Likelihoods are *not* measures, but density functions. They are callable
and also support the`DensityInterface` API. If `::AbstractLikelihood`, then

```julia
DensityInterface.DensityKind(ℒ) == IsDensity()
log(ℒ(θ)) ≈ logdensityof(ℒ, θ)
```

Given a transition kernel `k(θ)` (a function that takes a parameter object
and returns a measure) and an observation `x`, the recommended way to create a
likelihood object is

```julia
ℒ = likelihoodof(k, x)
ℒ isa AbstractLikelihood
```

Then

```julia
`log(ℒ(θ))` ≈ logdensityof(ℒ, θ) ≈ logdensityof(k(θ), x)
```

See [`likelihoodof`](@ref) for details on the mathematical semantics.
`k` and `x`.

To make likelihood-like types that are not subtypes of `AbstractLikelihood`
compatible with the `MeasureBase` likelihoods and Lebesque integrals,
specialize

```
MeasureBase.as_likelihood(ℒ::MyLikelihoodType) = likelihoodof(..., ...)
MeasureBase.as_integrand(ℒ::MyLikelihoodType) = AbstractLikelihood(ℒ)
```

Likelihood-like types that are not subtyles of `AbstractLikelihood` can made
compatible by specializing [`as_likelihood(L::MyLikelihoodType)`](@ref) and
[`as_integrand(L::MyLikelihoodType)`](@ref). By default, this is implemented
for objects like

```julia
L = Base.Fix2(logdensityof, x) ∘ f
L = FuncDensity(Base.Fix2(logdensityof, x) ∘ f)
L = LogFuncDensity(Base.Fix2(logdensityof, x) ∘ f)
```
"""
abstract type AbstractLikelihood <: Function end
export AbstractLikelihood

@inline AbstractLikelihood(l) = as_likelihood(L)::AbstractLikelihood

Base.convert(::Type{AbstractLikelihood}, l::AbstractLikelihood) = l
Base.convert(::Type{AbstractLikelihood}, l) = AbstractLikelihood(l)


"""
    likelihood_kernel(ℒ::AbstractLikelihood)

Return the Markov kernel that is part of likelihood `ℒ`.

If `ℒ = likelihoodof(k, x)` then `likelihood_kernel(ℒ)` must return an
equivalent of `k` (typically but not necessarily `k` itself).
"""
function likelihood_kernel end
export likelihood_kernel


"""
    likelihood_obs(ℒ::AbstractLikelihood)

Return the observation that is part of likelihood `ℒ`.

If `ℒ = likelihoodof(k, x)` then `likelihood_obs(ℒ)` must return an
equivalent of `x` (typically but not necessarily `x` itself).
```
"""
function likelihood_obs end
export likelihood_obs


"""
    MeasureBase.as_likelihood(L)::AbstractLikelihood

Turn a likelihood_like object `L` into an `AbstractLikelihood`.

Likelihood-like types that are not subtyles of `AbstractLikelihood` can made
compatible by specializing [`as_likelihood(L::MyLikelihoodType)`](@ref) and
[`as_integrand(L::MyLikelihoodType)`](@ref):

```julia
MeasureBase.as_likelihood(l::_SimpleLikelihood1) = likelihood_of(..., ...)
MeasureBase.as_integrand(l::_SimpleLikelihood1) = MeasureBase.as_likelihood(l)
```

By default, this is implemented for objects like

```julia
L = Base.Fix2(logdensityof, x) ∘ f
L = FuncDensity(Base.Fix2(logdensityof, x) ∘ f)
L = LogFuncDensity(Base.Fix2(logdensityof, x) ∘ f)
"""
function as_likelihood end
export as_likelihood

@inline as_likelihood(l::AbstractLikelihood) = l

@inline as_integrand(l::AbstractLikelihood) = l


(L::AbstractLikelihood)(p) = densityof(L, p)


DensityInterface.DensityKind(::AbstractLikelihood) = IsDensity()


_eval_k(L::AbstractLikelihood, p) = asmeasure(likelihood_kernel(L)(p))

function DensityInterface.logdensityof(L::AbstractLikelihood, p)
    logdensityof(_eval_k(L, p), likelihood_obs(L))
end

function DensityInterface.densityof(L::AbstractLikelihood, p)
    exp(ULogarithmic, logdensityof(_eval_k(L, p), likelihood_obs(L)))
end


const _SimpleLikelihood1 = ComposedFunction{<:Base.Fix2{typeof(densityof),<:Any},<:Any}
as_likelihood(l::_SimpleLikelihood1) = likelihood_of(l.inner, l.outer.x)
as_integrand(l::_SimpleLikelihood1) = as_likelihood(l)

const _SimpleLikelihood2 = DensityInterface.FuncDensity{<:ComposedFunction{<:Base.Fix2{typeof(densityof),<:Any},<:Any}}
as_likelihood(l::_SimpleLikelihood2) =  likelihood_of(l._f.inner, l._f.outer.x)
as_integrand(l::_SimpleLikelihood2) = as_likelihood(l)

const _SimpleLikelihood3 = DensityInterface.LogFuncDensity{<:ComposedFunction{<:Base.Fix2{typeof(logdensityof),<:Any},<:Any}}
as_likelihood(l::_SimpleLikelihood3) =  likelihood_of(l._log_f.inner, l._log_f.outer.x)
as_integrand(l::_SimpleLikelihood3) = as_likelihood(l)

const _SimpleLogLikelihood1 = ComposedFunction{<:Base.Fix2{typeof(logdensityof),<:Any},<:Any}
as_integrand_exp(l::_SimpleLogLikelihood1) = likelihood_of(l.inner, l.outer.x)



@doc raw"""
    struct Likelihood <: AbstractLikelihood

Default result of [`likelihoodof(k, x)`](@ref).

See [`AbstractLikelihood`](@ref) and [`likelihoodof`](@ref) for details.
"""
struct Likelihood{K,X} <: AbstractLikelihood
    k::K
    x::X

    Likelihood{K,X}(k, x) where {K,X} = new{K,X}(k, x)
end
export Likelihood

# For type stability, in case k is a type (resp. a constructor):
Likelihood(k, x::X) where {X} = Likelihood{Core.Typeof(k),X}(k, x)

likelihood_kernel(L::Likelihood) = L.k
likelihood_obs(L::Likelihood) = L.x

function Pretty.quoteof(L::Likelihood)
    k = Pretty.quoteof(L.k)
    x = Pretty.quoteof(L.x)
    :(Likelihood($k, $x))
end

function Base.show(io::IO, L::Likelihood)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, L)
end


# basemeasure(L::Likelihood) = @error "Likelihood requires local base measure"


@doc raw"""
    likelihoodof(k, x)::AbstractLikelihood

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

`ℒ = likelihoodof(k, x)` has the mathematical interpretation

```math
\mathcal{L}_x(\theta) = \frac{\rm{d}\, k(\theta)}{\rm{d}\, \chi}(x)
```

`likelihoodof` must return an object that implements the
[`DensityInterface`](https://github.com/JuliaMath/DensityInterface.jl)` API
and `ℒ = likelihoodof(k, x)` must satisfy

```julia
log(ℒ(θ)) == logdensityof(ℒ, θ) ≈ logdensityof(k(θ), x)

DensityKind(ℒ) isa IsDensity
```

[`likelihood_kernel(ℒ)`](@ref) must return an equivalent of `k` and
[`likelihood_obs(ℒ)`](@ref) must return an equivalent of `x` (typically, but
not necessarily, `k` and `x` themselves).

By default, an instance of [`MeasureBase.Likelihood`](@ref) is returned.
"""
function likelihoodof end
export likelihoodof

likelihoodof(k, x) = Likelihood(k, x)



###############################################################################
# At the least, we need to think through in some more detail whether
# (log-)likelihood ratios expressed in this way are correct and useful. For now
# this code is commented out; we may remove it entirely in the future.

# export log_likelihood_ratio

# """
#     log_likelihood_ratio(L::Likelihood, p, q)

# Compute the log of the likelihood ratio, in order to compare two choices for
# parameters. This is computed as

#     logdensity_rel(_eval_k(L, p), L.k(q), L.x)

# Since `logdensity_rel` can leave common base measure unevaluated, this can be
# more efficient than

#     logdensityof(_eval_k(L, p), L.x) - logdensityof(L.k(q), L.x)
# """
# log_likelihood_ratio(L::Likelihood, p, q) = logdensity_rel(_eval_k(L, p), L.k(q), L.x)

# # likelihoodof(k, x; kwargs...) = likelihoodof(k, x, NamedTuple(kwargs))

# export likelihood_ratio

# """
#     likelihood_ratio(L::Likelihood, p, q)

# Compute the log of the likelihood ratio, in order to compare two choices for
# parameters. This is equal to

#     density_rel(_eval_k(L, p), L.k(q), L.x)

# but is computed using LogarithmicNumbers.jl to avoid underflow and overflow.
# Since `density_rel` can leave common base measure unevaluated, this can be
# more efficient than

#     logdensityof(_eval_k(L, p), L.x) - logdensityof(L.k(q), L.x)
# """
# function likelihood_ratio(L::Likelihood, p, q)
#     exp(ULogarithmic, logdensity_rel(_eval_k(L, p), L.k(q), L.x))
# end
