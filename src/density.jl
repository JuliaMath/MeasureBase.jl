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

logdensityof(d::Density, x) = logdensity_rel(d.μ, d.base, x)

logdensity_def(d::Density, x) = logdensityof(d, x)

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

@inline function insupport(d::DensityMeasure, x)
    ifelse(insupport(d.base, x), logdensityof(d.f, x) > -Inf, false)
end

basemeasure(μ::DensityMeasure) = μ.base

logdensity_def(μ::DensityMeasure, x) = logdensityof(μ.f, x)

density_def(μ::DensityMeasure, x) = densityof(μ.f, x)


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

@inline function logdensityof(μ, x)
    t() = dynamic(unsafe_logdensityof(μ, x))
    f() = -Inf
    ifelse(insupport(μ, x), t, f)()
end

export unsafe_logdensityof

# https://discourse.julialang.org/t/counting-iterations-to-a-type-fixpoint/75876/10?u=cscherrer
@inline function unsafe_logdensityof(μ::M, x) where {M}
    ℓ_0 = logdensity_def(μ, x)
    b_0 = μ
    Base.Cartesian.@nexprs 10 i -> begin  # 10 is just some "big enough" number
        b_{i} = basemeasure(b_{i-1})
        if b_{i} isa typeof(b_{i-1})
            return ℓ_{i-1}
        end
        ℓ_{i} = let Δℓ_{i} = logdensity_def(b_{i}, x)
            ℓ_{i-1} + Δℓ_{i}
        end
    end
    return ℓ_10
end

export density_rel

@inline density_rel(μ, ν, x) = exp(logdensity_rel(μ, ν, x))

@inline function logdensity_rel(μ::M, ν::N, x::X) where {M,N,X}
    insupport(μ, x) || begin
        insupport(ν, x) || return NaN
        return -Inf
    end
    insupport(ν, x) || return Inf

    μ_depth = basemeasure_depth(μ)
    ν_depth = basemeasure_depth(ν)

    Δdepth = μ_depth - ν_depth

    if Δdepth > 0
        (ℓ₊, μ_0) = logdensity_steps(μ, x, Δdepth)
        ℓ₋ = zero(ℓ₊)
        ν_0 = ν
    elseif Δdepth < 0
        (ℓ₋, ν_0) = logdensity_steps(ν, x, -Δdepth)
        ℓ₊ = zero(ℓ₋)
        μ_0 = μ
    else
        ℓ₊ = ℓ₋ = 0.0
        μ_0 = μ
        ν_0 = ν
    end

    @assert basemeasure_depth(μ_0) == basemeasure_depth(ν_0)

    Base.Cartesian.@nexprs 10 i -> begin  # 10 is just some "big enough" number
        if μ_{i-1} == ν_{i-1}
            return ℓ₊ - ℓ₋
        elseif static_hasmethod(logdensity_def, Tuple{typeof(μ_{i-1}), typeof(ν_{i-1}), typeof(x)})
            return ℓ₊ - ℓ₋ + logdensity_def(μ_{i-1}, ν_{i-1}, x)
        end

        ℓ₊ += logdensity_def(μ_{i-1}, x)
        ℓ₋ += logdensity_def(ν_{i-1}, x)
        μ_{i} = basemeasure(μ_{i-1})
        ν_{i} = basemeasure(ν_{i-1})
    end
       
    @warn """
    No common base measure for
        $μ
    and
        $ν

    Returning a relative log-density of NaN. If this is incorrect, add a
    three-argument method
        logdensity_def(μ, ν, x)
    """
    return NaN
end

logdensity_steps(μ, x, ::StaticInt{0}) = (zero(logdensity_def(μ, x)), μ)

@inline function logdensity_steps(μ, x, ::StaticInt{n}) where {n}
    ℓ_0 = logdensity_def(μ, x)
    b_0 = μ
    
    Base.Cartesian.@nexprs 10 i -> begin  # 10 is just some "big enough" number
        i > n && begin
            return (ℓ_{i-1}, b_{i-1})
        end
        b_{i} = basemeasure(b_{i-1})
        ℓ_{i} = let Δℓ_{i} = logdensity_def(b_{i}, x)
            ℓ_{i-1} + Δℓ_{i}
        end       
    end
    return (ℓ_10, b_10)
end

@inline function _logdensity_rel(α::A, β::B, x::X, ℓ) where {A,B,X}
    if static_hasmethod(logdensity_def, Tuple{A,B,X})
        return ℓ + logdensity_def(α, β, x)
    elseif static_hasmethod(logdensity_def, Tuple{B,A,X})
        return ℓ + logdensity_def(β, α, x)
    else
        @warn """
        No method 
        logdensity(::$A, ::$B, ::$X)
        """
        return oftype(ℓ, NaN)
    end
end

export densityof
export logdensityof

export density_def

density_def(μ, ν::AbstractMeasure, x) = exp(logdensity_def(μ, ν, x))
density_def(μ, x) = exp(logdensity_def(μ, x))

"""
    rebase(μ, ν)

Express `μ` in terms of a density over `ν`. Satisfies
```
basemeasure(rebase(μ, ν)) == ν
density(rebase(μ, ν)) == 𝒹(μ,ν)
``` 
"""
rebase(μ, ν) = ∫(𝒹(μ,ν), ν)