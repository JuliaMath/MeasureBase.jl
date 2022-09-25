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
# Densities

"""
    struct Density{M,B}
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

"""
    struct DensityMeasure{F,B,L} <: AbstractMeasure
        density :: F
        base    :: B
    end

A `DensityMeasure` is a measure defined by a density with respect to some other
"base" measure 
"""
struct DensityMeasure{F,B,L} <: AbstractDensityMeasure
    f::F
    base::B
end


function Pretty.tile(μ::DensityMeasure{F,B}) where {F,B}
    result = Pretty.literal("DensityMeasure ∫(")
    result *= Pretty.pair_layout(Pretty.tile(μ.f), Pretty.tile(μ.base); sep = ", ")
    result *= Pretty.literal(")")
end


Base.log(d::Density) = logdensity(d.μ, d.base)

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


export 𝒹

"""
    𝒹(μ, base)

Compute the density (Radom-Nikodym derivative) of μ with respect to `base`.
"""
function 𝒹(μ, base)
    return density(μ, base)
end

density(μ, base) = Density(μ, base)

####################################################################################
# Log-densities

# TODO: Add methods for `exp ∘ (d::Density)` and `log ∘ (d::Density)`

"""
    struct LogDensity{M,B}
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

Base.exp(d::LogDensity) = density(d.μ, d.base)
Base.:∘(::typeof(exp), d::LogDensity) = density(d.μ, d.base)

logdensity(μ, base) = LogDensity(μ, base)

export log𝒹

"""
    log𝒹(μ, base)

Compute the log-density of μ with respect to `base`.
"""
log𝒹(μ, base) = logdensity(μ, base)

(f::LogDensity)(x) = logdensity_rel(f.μ, f.base, x)

(f::Density)(x) = density_rel(f.μ, f.base, x)

logdensity_def(d::Density, x) = logdensityof(d, x)

"""
    struct LogDensityMeasure{F,B,L} <: AbstractDensityMeasure
        density :: F
        base    :: B
    end

A `DensityMeasure` is a measure defined by a density with respect to some other
"base" measure 
"""
struct LogDensityMeasure{F,B} <: AbstractDensityMeasure
    f::F
    base::B
end

function Pretty.tile(μ::LogDensityMeasure{F,B}) where {F,B}
    result = Pretty.literal("LogDensityMeasure ∫exp(")
    result *= Pretty.pair_layout(Pretty.tile(μ.f), Pretty.tile(μ.base); sep = ", ")
    result *= Pretty.literal(")")
end

densitymeasure(f, base) = _densitymeasure(f, base, DensityKind(f))

_densitymeasure(f, base, ::IsDensity) = DensityMeasure(f, base)

function _densitymeasure(f, base, _)
    @error """
    The first argument of `DensityMeasure`" must be `::IsDensity`. To pass a
    function, first wrap it in `DensityInterface.funcdensity` or
    `DensityInterface.logfuncdensity`. 
    """
end

basemeasure(μ::DensityMeasure) = μ.base

logdensity_def(μ::DensityMeasure, x) = logdensityof(μ.f, x)

density_def(μ::DensityMeasure, x) = densityof(μ.f, x)

export ∫exp

"""
    ∫exp(f, base::AbstractMeasure)

Define a new measure in terms of a log-density `f` over some measure `base`.
"""
∫exp(f, base) = LogDensityMeasure(logfuncdensity(f), base)

∫exp(f::LogFuncDensity, base) = LogDensityMeasure(f, base)

function ∫exp(f::FuncDensity, base)
    @error "Can't call `∫exp` on a `FuncDensity`; use `∫` instead"
end
