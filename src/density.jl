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
struct Density{M,B,L}
    μ::M
    base::B
    log::L
end

export 𝒹

export log𝒹

log𝒹(μ, base) = Density(μ, base, Val{true}())

"""
    𝒹(μ::AbstractMeasure, base::AbstractMeasure; log=false)

Compute the Radom-Nikodym derivative (or its log, if `log=false`) of μ with
respect to `base`.
"""
function 𝒹(μ::AbstractMeasure, base::AbstractMeasure; log=false)
    return Density(μ, base, Val(log))
end

(f::Density{M,B,Val{true}})(x) where {M,B} = logdensity(f.μ, f.base, x)


(f::Density{M,B,Val{false}})(x) where {M,B} = density(f.μ, f.base, x)

"""
    struct DensityMeasure{F,B} <: AbstractMeasure
        density :: F
        base    :: B
    end

A `DensityMeasure` is a measure defined by a density with respect to some other
"base" measure 
"""
struct DensityMeasure{F,B,L} <: AbstractMeasure
    f    :: F
    base :: B
    log  :: L
end

function Base.show(io::IO, μ::DensityMeasure{F,B,Val{L}}) where {F,B,L}
    io = IOContext(io, :compact => true)
    print(io, "DensityMeasure ")
    print(io, "∫(", μ.f)
    print(io, ", ", μ.base)
    print(io, "; log = ", L, ")")
end

function Base.rand(rng::AbstractRNG, T::Type, d::DensityMeasure)
    x = rand(rng, T, d.base)
    WeightedMeasure(d.f(x), Dirac(x))
end

basemeasure(μ::DensityMeasure) = μ.base

logdensity(μ::DensityMeasure{F,B,Val{true}}, x) where {F,B} = μ.f(x)

density(μ::DensityMeasure{F,B,Val{false}}, x) where {F,B} = μ.f(x)

logdensity(μ::DensityMeasure{F,B,Val{false}}, x) where {F,B} = log(density(μ,x))

export ∫

"""
    ∫(f, base::AbstractMeasure)

Define a new measure in terms of a density `f` over some measure `base`.
"""
∫(f, base::AbstractMeasure) = DensityMeasure(f, base, Val(false))

∫(μ::AbstractMeasure, base::AbstractMeasure) = ∫exp(log𝒹(μ, base), base)


export ∫exp

"""
    ∫exp(f, base::AbstractMeasure; log=false)

Define a new measure in terms of a density `f` over some measure `base`.
"""
∫exp(f,μ) = DensityMeasure(f,μ,Val{true}())

# TODO: `density` and `logdensity` functions for `DensityMeasure`

function logdensity(μ::T, ν::T, x) where {T<:AbstractMeasure}
    μ==ν && return 0.0
    invoke(logdensity, Tuple{AbstractMeasure, AbstractMeasure, typeof(x)}, μ, ν, x)
end

function logdensity(μ::AbstractMeasure, ν::AbstractMeasure, x)
    α = basemeasure(μ)
    β = basemeasure(ν)

    # If α===μ and β===ν, The recursive call would be exactly the same as the
    # original one. We need to break the recursion.
    if α===μ && β===ν
        @warn """
        No method found for logdensity(μ, ν, x) where
        typeof(μ) == $(typeof(μ))
        typeof(ν) == $(typeof(ν))

        Returning NaN. If this is incorrect, please add a method        
        logdensity(μ::$(typeof(μ)), ν::$(typeof(ν)), x)
        """
        return NaN
    end

    # Infinite or NaN results occur when outside the support of α or β, 
    # and also when one measure is singular wrt the other. Computing the base
    # measures first is often much cheaper, and allows the numerically-intensive
    # computation to "fall through" in these cases.
    # TODO: Add tests to check that NaN cases work properly
    ℓ = logdensity(α, β, x)
    isfinite(ℓ) || return ℓ

    ℓ += logdensity(μ, x)
    ℓ -= logdensity(ν, x)

    return ℓ
end

function logpdf(d::AbstractMeasure, x)
    _logpdf(d, basemeasure(d), x, zero(Float64))
end

@inline function _logpdf(d::AbstractMeasure, β::AbstractMeasure, x, ℓ::Float64)
    d === β && return ℓ
    Δℓ = logdensity(d, x)
    # @show Δℓ, d
    _logpdf(β, basemeasure(β), x, ℓ + Δℓ)
end

logdensity(::Lebesgue, ::Lebesgue, x) = 0.0

# logdensity(::Lebesgue{ℝ}, ::Lebesgue{ℝ}, x) = zero(x)

export density

density(μ, ν::AbstractMeasure, x) = exp(logdensity(μ, ν, x))

density(μ, x) = exp(logdensity(μ, x))
