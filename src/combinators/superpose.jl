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
    components = (μ, ν)
    superpose(components)
end

logdensity_def(μ::SuperpositionMeasure, x) =
    logsumexp((logdensity_def(m, x) for m in μ.components))

basemeasure(μ::SuperpositionMeasure) = superpose(map(basemeasure, μ.components))

tbasemeasure_type(::Type{S}) where {C,S<:SuperpositionMeasure{C}}
    SuperpositionMeasure{tmap(tbasemeasure_type, C)}
end

# TODO: Fix `rand` method (this one is wrong)
# function Base.rand(μ::SuperpositionMeasure{X,N}) where {X,N}
#     return rand(rand(μ.components))
# end
