export SuperpositionMeasure

@doc raw"""
    struct SuperpositionMeasure{X,NT} <: AbstractMeasure
        components :: NT
    end
Superposition of measures is analogous to mixture distributions, but (because
measures need not be normalized) requires no scaling.
The superposition of two measures μ and ν can be more concisely written as μ + ν.
Superposition measures satisfy
    
    basemeasure(μ + ν) == basemeasure(μ) + basemeasure(ν)


```math
    \begin{aligned}\frac{\mathrm{d}(\mu+\nu)}{\mathrm{d}(\alpha+\beta)} & =\frac{f\,\mathrm{d}\alpha+g\,\mathrm{d}\beta}{\mathrm{d}\alpha+\mathrm{d}\beta}\\
     & =\frac{f\,\mathrm{d}\alpha}{\mathrm{d}\alpha+\mathrm{d}\beta}+\frac{g\,\mathrm{d}\beta}{\mathrm{d}\alpha+\mathrm{d}\beta}\\
     & =\frac{f}{1+\frac{\mathrm{d}\beta}{\mathrm{d}\alpha}}+\frac{g}{\frac{\mathrm{d}\alpha}{\mathrm{d}\beta}+1}\\
     & =\frac{f}{1+\left(\frac{\mathrm{d}\alpha}{\mathrm{d}\beta}\right)^{-1}}+\frac{g}{\frac{\mathrm{d}\alpha}{\mathrm{d}\beta}+1}\ .
    \end{aligned}
```
"""
struct SuperpositionMeasure{C} <: AbstractMeasure
    components::C
end

function Pretty.tile(d::SuperpositionMeasure)
    result = Pretty.literal("SuperpositionMeasure(")
    result *= Pretty.list_layout([Pretty.tile.(d.components)...])
    result *= Pretty.literal(")")
end

testvalue(μ::SuperpositionMeasure) = testvalue(first(μ.components))

# SuperpositionMeasure(ms :: AbstractMeasure...) = SuperpositionMeasure{X,length(ms)}(ms)

# SuperpositionMeasure(m::NTuple{N, Measure{X}}) where {N,X} = SuperpositionMeasure(m...)

# Base.length(::SuperpositionMeasure{X,N}) where {X,N} = N

# function Base.:+(μ::SuperpositionMeasure{N1}, ν::SuperpositionMeasure{N2}) where {N1,N2}
#     components = (μ.components..., ν.components...)
#     SuperpositionMeasure{X, N1+N2}(components)
# end

# function Base.:+(μ::AbstractMeasure, ν::SuperpositionMeasure{X,N}) where {X,N}
#     components = (μ, ν.components...)
#     SuperpositionMeasure{X,N+1}(components)
# end

# function Base.:+(μ::SuperpositionMeasure{X,N}, ν::AbstractMeasure) where {X,N}
#     components = (μ.components..., ν)
#     SuperpositionMeasure{X,N+1}(components)
# end

function Base.:+(μ::AbstractMeasure, ν::AbstractMeasure)
    superpose(μ, ν)
end

using LogarithmicNumbers

oneplus(x::ULogarithmic) = exp(ULogarithmic, log1pexp(x.log))

@inline function density_def(s::SuperpositionMeasure{Tuple{A,B}}, x) where {A,B}
    (μ, ν) = s.components

    # Jumping through some hoops here to avoid having this be a breaking
    # version. After the next breaking version, this can probably be changed to
    # the more familiar
    #     inμ || return exp(ULogarithmic, logdensity_def(ν, x))
    # etc
    let inμ = insupport(μ, x)
        if inμ isa False || !inμ
            return exp(ULogarithmic, logdensity_def(ν, x))
        end
    end
    let inν = insupport(ν, x)
        if !inν || inν isa False
            return exp(ULogarithmic, logdensity_def(μ, x))
        end
    end
    α = basemeasure(μ)
    β = basemeasure(ν)
    dμ_dα = exp(ULogarithmic, logdensity_def(μ, x))
    dν_dβ = exp(ULogarithmic, logdensity_def(ν, x))
    dα_dβ = exp(ULogarithmic, logdensity_rel(α, β, x))
    dβ_dα = inv(dα_dβ)
    return dμ_dα / oneplus(dβ_dα) + dν_dβ / oneplus(dα_dβ)
end

using LogExpFunctions

@inline function logdensity_def(μ::T, ν::T, x::Any) where T<:(SuperpositionMeasure{Tuple{A, B}} where {A, B})
    if μ === ν
        return zero(return_type(logdensity_def, (μ, x)))
    else
        return logdensity_def(μ,x) - logdensity_def(ν, x)
    end
end

@inline function logdensity_def(s::SuperpositionMeasure{Tuple{A,B}}, β, x) where {A,B}
    (μ, ν) = s.components
    insupport(μ, x) || return logdensity_rel(ν, β, x)
    insupport(ν, x) || return logdensity_rel(μ, β, x)
    return logaddexp(logdensity_rel(μ, β, x), logdensity_rel(ν, β, x))
end

@inline function logdensity_def(s::SuperpositionMeasure{Tuple{A,B}}, β::SuperpositionMeasure, x) where {A,B}
    (μ, ν) = s.components
    insupport(μ, x) || return logdensity_rel(ν, β, x)
    insupport(ν, x) || return logdensity_rel(μ, β, x)
    return logaddexp(logdensity_rel(μ, β, x), logdensity_rel(ν, β, x))
end

@inline function logdensity_def(s, β::SuperpositionMeasure{Tuple{A,B}}, x) where {A,B}
    -logdensity_def(β, s, x)
end

@inline logdensity_def(s::SuperpositionMeasure, x) = log(density_def(s, x))

basemeasure(μ::SuperpositionMeasure) = superpose(map(basemeasure, μ.components)...)

# TODO: Fix `rand` method (this one is wrong)
# function Base.rand(μ::SuperpositionMeasure{X,N}) where {X,N}
#     return rand(rand(μ.components))
# end

@inline function insupport(d::SuperpositionMeasure, x)
    any(d.components) do c
        dynamic(insupport(c, x))
    end
end