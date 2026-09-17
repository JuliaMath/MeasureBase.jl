
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

# Masks components outside of their support with -Inf:
@inline _masked_logd(ℓ, ins) = ifelse(_insupport_mask(ins), ℓ, oftype(ℓ, -Inf))

# Branch-free logsumexp over the components, valid for infinite entries:
@inline function _logsumexp_components(ℓs)
    m = reduce(max, ℓs)
    m_finite = ifelse(isfinite(m), m, zero(m))
    m_finite + log(sum(map(ℓ -> exp(ℓ - m_finite), ℓs)))
end

# The density of a superposition relative to the superposition of the
# component base measures, in log space: each component contributes its
# own density, divided by the density of the superposed base measures
# relative to its own base measure.
function logdensity_def(s::SuperpositionMeasure, x)
    cs = values(s.components)
    αs = map(basemeasure, cs)
    terms = map(cs, αs) do cᵢ, αᵢ
        ℓᵢ = _dynamic_logd(logdensity_def(cᵢ, x), x)
        log_dΣα_dαᵢ = _logsumexp_components(map(cs, αs) do cⱼ, αⱼ
            _masked_logd(logdensity_rel(αⱼ, αᵢ, x), insupport(cⱼ, x))
        end)
        _masked_logd(ℓᵢ - log_dΣα_dαᵢ, insupport(cᵢ, x))
    end
    _logsumexp_components(terms)
end

@inline function logdensity_rel_def(μ::T, ν::T, x) where {T<:SuperpositionMeasure}
    ℓ = logdensity_def(μ, x) - logdensity_def(ν, x)
    ifelse(μ === ν, zero(ℓ), ℓ)
end

function _superpos_logdensity_rel(s::SuperpositionMeasure, β, x)
    cs = values(s.components)
    ds = map(cs) do μ
        _masked_logd(logdensity_rel(μ, β, x), insupport(μ, x))
    end
    _logsumexp_components(ds)
end

@inline logdensity_rel_def(s::SuperpositionMeasure, β, x) = _superpos_logdensity_rel(s, β, x)

@inline logdensity_rel_def(s::SuperpositionMeasure, β::SuperpositionMeasure, x) =
    _superpos_logdensity_rel(s, β, x)

@inline logdensity_rel_def(s, β::SuperpositionMeasure, x) = -_superpos_logdensity_rel(β, s, x)

@inline density_def(s::SuperpositionMeasure, x) = exp(logdensity_def(s, x))

function basemeasure(μ::SuperpositionMeasure{<:Tuple})
    superpose(map(basemeasure, μ.components)...)
end

function basemeasure(μ::SuperpositionMeasure{<:AbstractArray})
    bases = map(basemeasure, μ.components)
    allequal(bases) ? weightedmeasure(log(length(bases)), first(bases)) : superpose(bases)
end

basemeasure(μ::SuperpositionMeasure) = superpose(map(basemeasure, μ.components))

function _component_masses(μ::SuperpositionMeasure)
    masses = map(massof, values(μ.components))
    total = sum(masses)
    total isa AbstractUnknownMass && throw(
        ArgumentError("Cannot sample from a superposition of measures of unknown mass"),
    )
    return map(dynamic, masses), dynamic(total)
end

function rand_impl(ctx::GenContext, μ::SuperpositionMeasure)
    components = values(μ.components)
    masses, total = _component_masses(μ)
    threshold = rand(get_rng(ctx), get_precision(ctx)) * total
    csum = zero(threshold)
    for (mass, c) in zip(masses, components)
        csum += mass
        csum >= threshold && return rand_impl(ctx, c)
    end
    return rand_impl(ctx, last(components))
end

# Batches of superpositions draw a batch from each component and select
# by mass, branch-free:
function batched_rand_impl(ctx::GenContext, μ::SuperpositionMeasure, sz::Dims)
    components = values(μ.components)
    masses, total = _component_masses(μ)
    thresholds = _rand_bulk(ctx, sz) .* total
    X = batched_rand_impl(ctx, first(components), sz)
    csum = first(masses)
    for (mass, c) in Iterators.drop(zip(masses, components), 1)
        X = ifelse.(thresholds .<= csum, X, batched_rand_impl(ctx, c, sz))
        csum += mass
    end
    return X
end

@inline function insupport(d::SuperpositionMeasure, x)
    mapreduce(c -> _insupport_mask(insupport(c, x)), |, values(d.components))
end


@inline mspace_flatsize(μ::SuperpositionMeasure) = mspace_flatsize(typeof(μ))
@inline mspace_flatsize(::Type{<:SuperpositionMeasure{C}}) where {C<:AbstractArray} = _scalar_or_unknown(mspace_flatsize(eltype(C)))
@inline mspace_flatsize(::Type{<:SuperpositionMeasure{C}}) where {C<:Tuple} = _common_scalar_flatsize(C)
@generated function _common_scalar_flatsize(::Type{C}) where {C<:Tuple}
    args = [:(mspace_flatsize($T)) for T in C.parameters]
    :(_all_scalar_sizes($(args...)))
end
@inline _all_scalar_sizes(::Tuple{}...) = ()
@inline _all_scalar_sizes(szs...) = NoMSpaceElementSize{typeof(szs)}()
