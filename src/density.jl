abstract type AbstractDensity end

@inline DensityKind(::AbstractDensity) = IsDensity()

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

export 𝒹

"""
    𝒹(μ::AbstractMeasure, base::AbstractMeasure)

Compute the Radom-Nikodym derivative of μ with respect to `base`.
"""
function 𝒹(μ::AbstractMeasure, base::AbstractMeasure)
    return Density(μ, base)
end

logdensityof(d::Density, x) = logdensityof(d.μ, d.base, x)

logdensityof(d::Density, x) = exp(logdensityof(d.μ, d.base, x))


"""
    struct DensityMeasure{F,B} <: AbstractMeasure
        density :: F
        base    :: B
    end

A `DensityMeasure` is a measure defined by a density with respect to some other
"base" measure 
"""
struct DensityMeasure{F,B} <: AbstractMeasure
    f::F
    base::B
end

function Pretty.tile(μ::DensityMeasure{F,B}) where {F,B}
    result = Pretty.literal("DensityMeasure ∫(")
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

basemeasure_type(::Type{DensityMeasure{F,B}}) where {F,B} = B

logdensity_def(μ::DensityMeasure, x) = logdensityof(μ.f, x)

density_def(μ::DensityMeasure, x) = densityof(μ.f, x)

densityof(μ::AbstractMeasure, ν::AbstractMeasure, x) = exp(logdensityof(μ, ν, x))

export ∫

"""
    ∫(f, base::AbstractMeasure)

Define a new measure in terms of a density `f` over some measure `base`.
"""
∫(f::Function, base::AbstractMeasure) = DensityMeasure(funcdensity(f), base)

∫(f, base::AbstractMeasure) = _densitymeasure(f, base, DensityKind(f))

# ∫(μ::AbstractMeasure, base::AbstractMeasure) = ∫(𝒹(μ, base), base)

export ∫exp

"""
    ∫exp(f, base::AbstractMeasure)

Define a new measure in terms of a log-density `f` over some measure `base`.
"""
∫exp(f::Function, μ) = ∫(logfuncdensity(f), μ)

# TODO: `density` and `logdensity` functions for `DensityMeasure`

@inline function logdensityof(μ::T, ν::T, x) where {T<:AbstractMeasure}
    μ == ν && return 0.0
    invoke(logdensityof, Tuple{AbstractMeasure,AbstractMeasure,typeof(x)}, μ, ν, x)
end

@inline function logdensityof(μ::AbstractMeasure, ν::AbstractMeasure, x)
    α = basemeasure(μ)
    β = basemeasure(ν)

    # If α===μ and β===ν, The recursive call would be exactly the same as the
    # original one. We need to break the recursion.
    if α === μ && β === ν
        @warn """
        No method found for logdensity_def(μ, ν, x) where
        typeof(μ) == $(typeof(μ))
        typeof(ν) == $(typeof(ν))

        Returning NaN. If this is incorrect, please add a method        
        logdensity_def(μ::$(typeof(μ)), ν::$(typeof(ν)), x)
        """
        return NaN
    end

    # Infinite or NaN results occur when outside the support of α or β, 
    # and also when one measure is singular wrt the other. Computing the base
    # measures first is often much cheaper, and allows the numerically-intensive
    # computation to "fall through" in these cases.
    # TODO: Add tests to check that NaN cases work properly
    ℓ = logdensityof(α, β, x)
    isnan(ℓ) && return ℓ

    ℓ += logdensity_def(μ, x)
    ℓ -= logdensity_def(ν, x)

    return ℓ
end

@inline logdensityof(μ, x) = first(_logdensityof(μ, x))

# Because it's sometimes useful, this returns a pair (ℓ,r) where
# • ℓ is the log-density
# • r is the rootmeasure of μ
@inline function _logdensityof(μ, x)
    n = basemeasure_depth(proxy(μ)) - static(1)
    
    β = basemeasure(μ, x)
    ℓ = logdensity_def(μ, x)
    return _logdensityof(β, x, ℓ, n)
end

@generated function _logdensityof(μ, x, ℓ, ::StaticInt{n}) where {n}
    nsteps = max(n, 0)
    quote
        $(Expr(:meta, :inline))
        Base.Cartesian.@nexprs $nsteps i -> begin
            Δℓ = logdensity_def(μ, x)
            μ = basemeasure(μ, x)
            ℓ += Δℓ
        end
        return (ℓ,μ)
    end
end

# logdensity_def(::Lebesgue{ℝ}, ::Lebesgue{ℝ}, x) = zero(x)

export densityof
export logdensityof

export density_def


density_def(μ, ν::AbstractMeasure, x) = exp(logdensity_def(μ, ν, x))

density_def(μ, x) = exp(logdensity_def(μ, x))
