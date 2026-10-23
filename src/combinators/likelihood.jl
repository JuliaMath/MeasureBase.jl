"""
    abstract type AbstractLikelihood <: Function

Abstract supertype for likelihood objects.

Likelihoods are *not* measures, but density functions. They are callable
and also support the DensityInterface API. If `ℒ isa AbstractLikelihood`,
then

```julia
DensityInterface.DensityKind(ℒ) == IsDensity()
log(ℒ(θ)) ≈ logdensityof(ℒ, θ)
```

Given a transition kernel `k(θ)` (a function that takes a parameter object
and returns a measure) and an observation `x`, the recommended way to create
a likelihood object is

```julia
ℒ = likelihoodof(k, x)
ℒ isa AbstractLikelihood
```

Then

```julia
log(ℒ(θ)) ≈ logdensityof(ℒ, θ) ≈ logdensityof(k(θ), x)
```

See [`likelihoodof`](@ref) for details on the mathematical semantics of
`k` and `x`.

Likelihood-like types that are not subtypes of `AbstractLikelihood` can be
made compatible with the `MeasureBase` likelihoods and Lebesgue integrals by
specializing [`MeasureBase.as_likelihood`](@ref) and
[`MeasureBase.as_integrand`](@ref).
"""
abstract type AbstractLikelihood <: Function end
export AbstractLikelihood

@inline AbstractLikelihood(l) = as_likelihood(l)::AbstractLikelihood

Base.convert(::Type{AbstractLikelihood}, l::AbstractLikelihood) = l
Base.convert(::Type{AbstractLikelihood}, l) = AbstractLikelihood(l)


"""
    likelihood_kernel(ℒ::AbstractLikelihood)

Return the transition kernel that is part of likelihood `ℒ`.

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
"""
function likelihood_obs end
export likelihood_obs


"""
    MeasureBase.as_likelihood(l)::AbstractLikelihood

Turn a likelihood-like object `l` into an `AbstractLikelihood`.

Likelihood-like types that are not subtypes of `AbstractLikelihood` can be
made compatible by specializing

```julia
MeasureBase.as_likelihood(l::MyLikelihoodType) = likelihoodof(..., ...)
MeasureBase.as_integrand(l::MyLikelihoodType) = MeasureBase.as_likelihood(l)
```

By default, this is implemented for objects like

```julia
l = Base.Fix2(densityof, x) ∘ f
l = FuncDensity(Base.Fix2(densityof, x) ∘ f)
l = LogFuncDensity(Base.Fix2(logdensityof, x) ∘ f)
```
"""
function as_likelihood end
export as_likelihood

@inline as_likelihood(l::AbstractLikelihood) = l

@inline as_integrand(l::AbstractLikelihood) = l


(ℒ::AbstractLikelihood)(p) = densityof(ℒ, p)


DensityInterface.DensityKind(::AbstractLikelihood) = IsDensity()


_eval_k(ℒ::AbstractLikelihood, p) = asmeasure(likelihood_kernel(ℒ)(p))

function DensityInterface.logdensityof(ℒ::AbstractLikelihood, p)
    logdensityof(_eval_k(ℒ, p), likelihood_obs(ℒ))
end

function DensityInterface.densityof(ℒ::AbstractLikelihood, p)
    exp(ULogarithmic, logdensityof(_eval_k(ℒ, p), likelihood_obs(ℒ)))
end


# Broadcasting over arrays of parameter values evaluates the measures the
# kernel returns for them as a parameter batch (see `logdensities` for
# arrays of measures), scalar observations in a single broadcast over the
# parameters (their columns for struct arrays):
const _LikelihoodArg = Union{AbstractLikelihood,Base.RefValue{<:AbstractLikelihood}}
@inline _unref(ℒ::AbstractLikelihood) = ℒ
@inline _unref(ℒ::Base.RefValue{<:AbstractLikelihood}) = ℒ[]

function Broadcast.broadcasted(::typeof(logdensityof), ℒ::_LikelihoodArg, Θ::AbstractArray)
    _likelihood_bc(_unref(ℒ), Θ, likelihood_obs(_unref(ℒ)))
end
function Broadcast.broadcasted(::typeof(densityof), ℒ::_LikelihoodArg, Θ::AbstractArray)
    Broadcast.broadcasted(Base.Fix1(exp, ULogarithmic), Broadcast.broadcasted(logdensityof, ℒ, Θ))
end

@inline _likelihood_bc(ℒ, Θ, x::Number) = _marginal_broadcast(_KernelLogd(likelihood_kernel(ℒ)), Θ, x)
@inline _likelihood_bc(ℒ, Θ, x) = Broadcast.broadcasted(logdensityof, _kernel_measures(ℒ, Θ), Ref(x))

struct _KernelMeasure{K} <: Function
    k::K
end
_KernelMeasure(k::Type{K}) where {K} = _KernelMeasure{Type{K}}(k)
@inline (f::_KernelMeasure)(θ) = asmeasure(f.k(θ))

struct _KernelLogd{K} <: Function
    k::K
end
_KernelLogd(k::Type{K}) where {K} = _KernelLogd{Type{K}}(k)
@inline (f::_KernelLogd)(θ, x) = _dynamic_logd(logdensityof(asmeasure(f.k(θ)), x), x)

# Broadcasting over struct arrays of parameters gives struct arrays of
# measures:
_kernel_measures(ℒ::AbstractLikelihood, Θ::AbstractArray) = broadcast(_KernelMeasure(likelihood_kernel(ℒ)), Θ)

# Integrands of density measures over arrays of variates:
_integrand_lds(ℒ::AbstractLikelihood, Θ::AbstractArray) = Broadcast.materialize(Broadcast.broadcasted(logdensityof, ℒ, Θ))


const _SimpleLikelihood1 = ComposedFunction{<:Base.Fix2{typeof(densityof),<:Any},<:Any}
as_likelihood(l::_SimpleLikelihood1) = likelihoodof(l.inner, l.outer.x)
as_integrand(l::_SimpleLikelihood1) = as_likelihood(l)

const _SimpleLikelihood2 = DensityInterface.FuncDensity{
    <:ComposedFunction{<:Base.Fix2{typeof(densityof),<:Any},<:Any},
}
as_likelihood(l::_SimpleLikelihood2) = likelihoodof(l._f.inner, l._f.outer.x)
as_integrand(l::_SimpleLikelihood2) = as_likelihood(l)

const _SimpleLikelihood3 = DensityInterface.LogFuncDensity{
    <:ComposedFunction{<:Base.Fix2{typeof(logdensityof),<:Any},<:Any},
}
as_likelihood(l::_SimpleLikelihood3) = likelihoodof(l._log_f.inner, l._log_f.outer.x)
as_integrand(l::_SimpleLikelihood3) = as_likelihood(l)

const _SimpleLogLikelihood1 = ComposedFunction{<:Base.Fix2{typeof(logdensityof),<:Any},<:Any}
as_integrand_exp(l::_SimpleLogLikelihood1) = likelihoodof(l.inner, l.outer.x)


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

likelihood_kernel(ℒ::Likelihood) = ℒ.k
likelihood_obs(ℒ::Likelihood) = ℒ.x

function Pretty.quoteof(ℒ::Likelihood)
    k = Pretty.quoteof(ℒ.k)
    x = Pretty.quoteof(ℒ.x)
    :(Likelihood($k, $x))
end

function Base.show(io::IO, ℒ::Likelihood)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, ℒ)
end


@doc raw"""
    likelihoodof(k, x)::AbstractLikelihood

Returns the likelihood of observing `x` under a family of probability
measures that is generated by a transition kernel `k(θ)`.

`k(θ)` maps points in the parameter space to measures (resp. objects that can
be converted to measures) on an implicit set `Χ` that contains values like
`x`.

`likelihoodof(k, x)` returns a likelihood object. A likelihood is **not** a
measure, it is a function from the parameter space to `ℝ₊`. Likelihood
objects can also be interpreted as "generic densities" (but **not** as
probability densities).

`likelihoodof(k, x)` implicitly chooses `ξ = rootmeasure(k(θ))` as the
reference measure on the observation set `Χ`. Note that this implicit
`ξ` **must** be independent of `θ`.

`ℒ = likelihoodof(k, x)` has the mathematical interpretation

```math
\mathcal{L}_x(\theta) = \frac{\rm{d}\, k(\theta)}{\rm{d}\, \chi}(x)
```

`likelihoodof` must return an object that implements the
[`DensityInterface`](https://github.com/JuliaMath/DensityInterface.jl) API
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
