"""
    abstract type MeasureBase.StdMeasure

Abstract supertype for standard measures.

Standard measures must be singleton types that represent common fundamental
measures such as [`StdUniform`](@ref), [`StdExponential`](@ref),
[`StdNormal`](@ref) and [`StdLogistic`](@ref).

A standard measure ([`StdUniform`](@ref), [`StdExponential`](@ref) and
[`StdNormal`](@ref)) is defined for every common Julia random number
generation function:

```
StdMeasure(rand) == StdUniform()
StdMeasure(randexp) == StdExponential()
StdMeasure(randn) == StdNormal()
```

[`StdLogistic`](@ref) has no associated random number generation function.

All standard measures must be normalized, i.e. [`massof`](@ref) always
returns one.
"""
abstract type StdMeasure <: AbstractMeasure end

@inline massof(::StdMeasure) = static(true)
@inline getdof(::StdMeasure) = static(1)

StdMeasure(::typeof(rand)) = StdUniform()
StdMeasure(::typeof(randexp)) = StdExponential()
StdMeasure(::typeof(randn)) = StdNormal()

@inline check_dof(::StdMeasure, ::StdMeasure) = nothing


# Transport between two equal standard measures:

@inline transport_def(::MU, μ::MU, x) where {MU<:StdMeasure} = x


# Transport between a standard measure and Dirac:

@inline transport_from_mvstd_with_rest(ν::Dirac, ::StdMeasure, x::Any) = ν.x, x

@inline transport_to_mvstd(ν::PowerMeasure{<:StdMeasure}, ::Dirac, ::Any) = Zeros{Bool}(map(_ -> 0, ν.axes))
