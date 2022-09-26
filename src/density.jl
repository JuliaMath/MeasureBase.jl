###################################################################
# Abstract types and methods

abstract type AbstractDensity <: Function end

@inline DensityKind(::AbstractDensity) = IsDensity()

abstract type AbstractDensityMeasure <: AbstractMeasure end

@inline function insupport(d::AbstractDensityMeasure, x)
    insupport(d.base, x) == true && isfinite(logdensityof(d.f, x))
end

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

Base.:∘(::typeof(log), d::Density) = logdensity(d.μ, d.base)

Base.log(d::Density) = log ∘ d

export 𝒹

"""
    𝒹(μ, base)

Compute the density (Radom-Nikodym derivative) of μ with respect to `base`.
"""
function 𝒹(μ, base)
    return density(μ, base)
end

density(μ, base) = Density(μ, base)

(f::Density)(x) = density_rel(f.μ, f.base, x)

#######################################################################################
# DensityMeasure

"""
    struct DensityMeasure{F,B,L} <: AbstractDensityMeasure
        density :: F
        base    :: B
    end

A `DensityMeasure` is a measure defined by a density with respect to some other
"base" measure 
"""
struct DensityMeasure{F,B} <: AbstractDensityMeasure
    f::F
    base::B
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
∫(f, base) = DensityMeasure(funcdensity(f), base)

∫(f::FuncDensity, base) = DensityMeasure(f, base)

function ∫(f::LogFuncDensity, base)
    @error "Can't call `∫` on a `LogFuncDensity`; use `∫exp` instead"
end

basemeasure(μ::DensityMeasure) = μ.base

logdensity_def(μ::DensityMeasure, x) = logdensityof(μ.f, x)

density_def(μ::DensityMeasure, x) = densityof(μ.f, x)

####################################################################################
# LogDensity

"""
    struct LogDensity{M,B} <: AbstractLogDensity
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
struct LogDensity{M,B} <: AbstractLogDensity
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
    return logdensity(μ, base)
end

logdensity(μ, base) = LogDensity(μ, base)

(f::LogDensity)(x) = logdensity_rel(f.μ, f.base, x)

#######################################################################################
# LogDensityMeasure

"""
    struct LogDensityMeasure{F,B,L} <: AbstractLogDensityMeasure
        density :: F
        base    :: B
    end

A `LogDensityMeasure` is a measure defined by a density with respect to some other
"base" measure 
"""
struct LogDensityMeasure{F,B} <: AbstractLogDensityMeasure
    f::F
    base::B
end

function Pretty.tile(μ::LogDensityMeasure{F,B}) where {F,B}
    result = Pretty.literal("LogDensityMeasure ∫(")
    result *= Pretty.pair_layout(Pretty.tile(μ.f), Pretty.tile(μ.base); sep = ", ")
    result *= Pretty.literal(")")
end

export ∫exp

"""
    ∫exp(f, base::AbstractMeasure)

Define a new measure in terms of a density `f` over some measure `base`.
"""
∫exp(f, base) = LogDensityMeasure(logfuncdensity(f), base)

∫exp(f::LogFuncDensity, base) = LogDensityMeasure(f, base)

function ∫exp(f::FuncDensity, base)
    @error "Can't call `∫exp` on a `FuncDensity`; use `∫` instead"
end

basemeasure(μ::LogDensityMeasure) = μ.base

logdensity_def(μ::LogDensityMeasure, x) = logdensityof(μ.f, x)

density_def(μ::LogDensityMeasure, x) = densityof(μ.f, x)
