
using LogarithmicNumbers
using LogExpFunctions

export SuperpositionMeasure

abstract type AbstractSuperpositionMeasure <: AbstractMeasure end

@doc raw"""
    struct SuperpositionMeasure{NT} <: AbstractMeasure
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
struct SuperpositionMeasure{C} <: AbstractSuperpositionMeasure
    components::C
end

massof(m::SuperpositionMeasure) = sum(massof, m.components)

function Pretty.tile(d::SuperpositionMeasure)
    result = Pretty.literal("SuperpositionMeasure(")
    result *= Pretty.list_layout([Pretty.tile.(d.components)...])
    result *= Pretty.literal(")")
end

testvalue(::Type{T}, μ::SuperpositionMeasure) where {T} = testvalue(T, first(μ.components))

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

@inline _ulogexp(x) = exp(ULogarithmic, dynamic(x))

function density_def(s::SuperpositionMeasure, x)
    cs = values(s.components)
    αs = map(basemeasure, cs)
    idxs = eachindex(cs)
    sum(idxs) do i
        dμᵢ_dαᵢ = _ulogexp(logdensity_def(cs[i], x))
        istrue(insupport(cs[i], x)) || return zero(dμᵢ_dαᵢ)
        dΣα_dαᵢ = sum(idxs) do j
            dαⱼ_dαᵢ = _ulogexp(logdensity_rel(αs[j], αs[i], x))
            istrue(insupport(cs[j], x)) ? dαⱼ_dαᵢ : zero(dαⱼ_dαᵢ)
        end
        dμᵢ_dαᵢ / dΣα_dαᵢ
    end
end

@inline function logdensity_def(μ::T, ν::T, x) where {T<:SuperpositionMeasure}
    if μ === ν
        return zero(return_type(logdensity_def, (μ, x)))
    else
        return logdensity_def(μ, x) - logdensity_def(ν, x)
    end
end

function _superpos_logdensity_rel(s::SuperpositionMeasure, β, x)
    cs = values(s.components)
    ds = map(cs) do μ
        istrue(insupport(μ, x)) ? dynamic(logdensity_rel(μ, β, x)) : -Inf
    end
    logsumexp(ds)
end

@inline logdensity_def(s::SuperpositionMeasure, β, x) = _superpos_logdensity_rel(s, β, x)

@inline logdensity_def(s::SuperpositionMeasure, β::SuperpositionMeasure, x) =
    _superpos_logdensity_rel(s, β, x)

@inline logdensity_def(s, β::SuperpositionMeasure, x) = -_superpos_logdensity_rel(β, s, x)

@inline logdensity_def(s::SuperpositionMeasure, x) = log(density_def(s, x))

function basemeasure(μ::SuperpositionMeasure{<:Tuple})
    superpose(map(basemeasure, μ.components)...)
end

function basemeasure(μ::SuperpositionMeasure{<:AbstractArray})
    bases = map(basemeasure, μ.components)
    allequal(bases) ? weightedmeasure(log(length(bases)), first(bases)) : superpose(bases)
end

basemeasure(μ::SuperpositionMeasure) = superpose(map(basemeasure, μ.components))

function Base.rand(rng::AbstractRNG, ::Type{T}, μ::SuperpositionMeasure) where {T}
    components = values(μ.components)
    masses = map(massof, components)
    total = sum(masses)
    total isa AbstractUnknownMass && throw(
        ArgumentError("Cannot sample from a superposition of measures of unknown mass"),
    )
    threshold = rand(rng) * dynamic(total)
    csum = zero(threshold)
    for (mass, c) in zip(masses, components)
        csum += dynamic(mass)
        csum >= threshold && return rand(rng, T, c)
    end
    return rand(rng, T, last(components))
end

@inline function insupport(d::SuperpositionMeasure, x)
    any(d.components) do c
        dynamic(insupport(c, x))
    end
end
