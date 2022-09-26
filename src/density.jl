###################################################################
# Abstract types and methods

abstract type AbstractDensity <: Function end

@inline DensityKind(::AbstractDensity) = IsDensity()

logdensityof(d::AbstractDensity, x) = logdensity_rel(d.μ, d.base, x)
densityof(d::AbstractDensity, x) = density_rel(d.μ, d.base, x)

####################################################################################
# Density

"""
    struct Density{M,B} <: AbstractDensity
        μ::M
        base::B
    end

For measures μ and ν with μ≪ν, the density of μ with respect to ν (also called
the Radon-Nikodym derivative dμ/dν) is a function f defined on the support of ν
with the property that for any measurable a ⊂ supp(ν), μ(a) = ∫ₐ f dν.
    
Because this function is often difficult to express in closed form, there are
many different ways of computing it. We therefore provide a formal
representation to allow comptuational flexibilty.
"""
struct Density{M,B} <: AbstractDensity
    μ::M
    base::B
end

Base.:∘(::typeof(log), d::Density) = logdensity_rel(d.μ, d.base)

Base.log(d::Density) = log ∘ d

export 𝒹

"""
    𝒹(μ, base)

Compute the density (Radom-Nikodym derivative) of μ with respect to `base`.
"""
function 𝒹(μ, base)
    return density_rel(μ, base)
end

density_rel(μ, base) = Density(μ, base)

(f::Density)(x) = density_rel(f.μ, f.base, x)

####################################################################################
# LogDensity

"""
    struct LogDensity{M,B} <: AbstractDensity
        μ::M
        base::B
    end

For measures μ and ν with μ≪ν, the density of μ with respect to ν (also called
the Radon-Nikodym derivative dμ/dν) is a function f defined on the support of ν
with the property that for any measurable a ⊂ supp(ν), μ(a) = ∫ₐ f dν.
    
Because this function is often difficult to express in closed form, there are
many different ways of computing it. We therefore provide a formal
representation to allow comptuational flexibilty.
"""
struct LogDensity{M,B} <: AbstractDensity
    μ::M
    base::B
end

Base.:∘(::typeof(exp), d::LogDensity) = density(d.μ, d.base)

Base.exp(d::LogDensity) = exp ∘ d

export log𝒹

"""
    log𝒹(μ, base)

Compute the density (Radom-Nikodym derivative) of μ with respect to `base`.
"""
function log𝒹(μ, base)
    return logdensity_rel(μ, base)
end

logdensity_rel(μ, base) = LogDensity(μ, base)

(f::LogDensity)(x) = logdensity_rel(f.μ, f.base, x)

#######################################################################################
# DensityMeasure

"""
    struct DensityMeasure{F,B} <: AbstractDensityMeasure
        density :: F
        base    :: B
    end

A `DensityMeasure` is a measure defined by a density or log-density with respect
to some other "base" measure.

Users should not call `DensityMeasure` directly, but should instead call `∫(f,
base)` (if `f` is a density) or `∫exp(f, base)` (if `f` is a log-density).
"""
struct DensityMeasure{F,B} <: AbstractMeasure
    f::F
    base::B

    function DensityMeasure(f::F, base::B) where {F,B}
        @assert DensityKind(f) isa IsDensity
        new{F,B}(f, base)
    end
end

function Base.propertynames(::DensityMeasure)
    (:density, :logdensity, :base)
end

function Base.getproperty(μ::DensityMeasure, s::Symbol)
    f = getfield(μ, :f)

    if s == :density
        return x -> densityof(f, x)
    elseif s == :logdensity
        return x -> logdensityof(f, x)
    elseif s == :base
        return getfield(μ, :base)
    end
end

@inline function insupport(d::DensityMeasure, x)
    insupport(d.base, x) == true && isfinite(logdensityof(d.f, x))
end

function Pretty.tile(μ::DensityMeasure{F,B}) where {F,B}
    result = Pretty.literal("DensityMeasure ∫(")
    result *= Pretty.pair_layout(Pretty.tile(μ.f), Pretty.tile(μ.base); sep = ", ")
    result *= Pretty.literal(")")
end

export ∫

"""
    ∫(f, base::AbstractMeasure)

Define a new measure in terms of a density `f` over some measure `base`.
"""
∫(f, base) = _densitymeasure(f, base, DensityKind(f))

_densitymeasure(f, base, ::IsDensity) = DensityMeasure(f, base)
_densitymeasure(f, base, ::HasDensity) = @error "..."
_densitymeasure(f, base, ::NoDensity) = DensityMeasure(funcdensity(f), base)

export ∫exp

"""
    ∫exp(f, base::AbstractMeasure)

Define a new measure in terms of a log-density `f` over some measure `base`.
"""
∫exp(f, base) = _logdensitymeasure(f, base, DensityKind(f))

_logdensitymeasure(f, base, ::IsDensity) = DensityMeasure(f, base)
_logdensitymeasure(f, base, ::HasDensity) = @error "..."
_logdensitymeasure(f, base, ::NoDensity) = DensityMeasure(logfuncdensity(f), base)

basemeasure(μ::DensityMeasure) = μ.base

logdensity_def(μ::DensityMeasure, x) = logdensityof(getfield(μ, :f), x)

density_def(μ::DensityMeasure, x) = densityof(getfield(μ, :f), x)
