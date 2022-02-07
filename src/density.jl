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

logdensityof(d::Density, x) = logdensityof(d.μ, x) - logdensityof(d.base, x)

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
    dynamic(insupport(μ, x)) || return -Inf
    return unsafe_logdensityof(μ, x)
end

export unsafe_logdensityof

# https://discourse.julialang.org/t/counting-iterations-to-a-type-fixpoint/75876/10?u=cscherrer
@inline function unsafe_logdensityof(μ::M, x) where {M}
    ℓ_0 = (logdensity_def(μ, x))
    b_0 = μ
    Base.Cartesian.@nexprs 10 i -> begin  # 10 is just some "big enough" number
        b_{i} = basemeasure(b_{i-1})
        b = b_{i}
        # @show b
        if b_{i} isa typeof(b_{i-1})
            return ℓ_{i-1}
        end
        ℓ_{i} = let Δt = @elapsed Δℓ = (logdensity_def(b_{i}, x))
            # @show Δt
            # @show Δℓ
            ℓ_{i-1} + Δℓ
        end
    end
    return ℓ_10
end

# # https://discourse.julialang.org/t/counting-iterations-to-a-type-fixpoint/75876/10?u=cscherrer
# @inline function unsafe_logdensityof(μ::M, x) where {M}
#     unsafe_logdensityof(μ, x, basemeasure_depth(μ))
# end

# @generated function unsafe_logdensityof(μ, x, ::StaticInt{N}) where {N}
#     q = quote
#         ℓ = logdensity_def(μ, x)
#     end

#     for j in 1:N
#         push!(q.args, quote
#             μ = basemeasure(μ)
#             Δℓ = logdensity_def(μ, x)
#             ℓ = ℓ + Δℓ
#         end)
#     end
#     return q
# end
   

@inline function logdensity_rel(μ::M, ν::N, x::X) where {M,N,X}
    (ℓ₊, α) = _logdensityof(μ, basemeasure(μ), x)
    (ℓ₋, β) = _logdensityof(ν, basemeasure(ν), x)
    return _logdensity_rel(α, β, x, ℓ₊ - ℓ₋)
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

# logdensity_def(::Lebesgue{ℝ}, ::Lebesgue{ℝ}, x) = zero(x)

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