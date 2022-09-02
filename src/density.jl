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
struct Density{M,B,L} <: AbstractDensity
    μ::M
    base::B
    log::L
end

Base.exp(d::Density{M,B,True})  = density(d.μ, d.base, False())
Base.log(d::Density{M,B,False}) = density(d.μ, d.base, True())

# TODO: Add methods for `exp ∘ (d::Density)` and `log ∘ (d::Density)`

density(μ, base, log) = Density(μ, base, log)

export 𝒹

export log𝒹

"""
    log𝒹(μ, base)

Compute the log-density of μ with respect to `base`.
"""
log𝒹(μ, base) = density(μ, base, True())

"""
    𝒹(μ, base)

Compute the density (Radom-Nikodym derivative) of μ with respect to `base`.
"""
function 𝒹(μ, base)
    return density(μ, base, False())
end

(f::Density{M,B,True})(x) where {M,B} = logdensity_rel(f.μ, f.base, x)

(f::Density{M,B,False})(x) where {M,B} = density_rel(f.μ, f.base, x)

logdensityof(d::Density, x) = logdensity_rel(d.μ, d.base, x)

logdensity_def(d::Density, x) = logdensityof(d, x)

"""
    struct DensityMeasure{F,B,L} <: AbstractMeasure
        density :: F
        base    :: B
        log     :: L
    end

A `DensityMeasure` is a measure defined by a density with respect to some other
"base" measure 
"""
struct DensityMeasure{F,B,L} <: AbstractMeasure
    f::F
    base::B
    log::L
end

function Pretty.tile(μ::DensityMeasure{F,B,True}) where {F,B}
    result = Pretty.literal("DensityMeasure ∫(")
    result *= Pretty.pair_layout(Pretty.tile(μ.f), Pretty.tile(μ.base); sep = ", ")
    result *= Pretty.literal("; log = true)")
end

function Pretty.tile(μ::DensityMeasure{F,B,False}) where {F,B}
    result = Pretty.literal("DensityMeasure ∫(")
    result *= Pretty.pair_layout(Pretty.tile(μ.f), Pretty.tile(μ.base); sep = ", ")
    result *= Pretty.literal("; log = false)")
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
    insupport(d.base, x) == true && isfinite(logdensityof(d.f, x))
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

∫exp(f, base::AbstractMeasure) = _densitymeasure(f, base, DensityKind(f))

"""
    logdensityof(m::AbstractMeasure, x) 

Compute the log-density of the measure `m` at `x`. Density is always relative,
but `DensityInterface.jl` does not account for this. For compatibility with
this, `logdensityof` for a measure is always implicitly relative to
[`rootmeasure(x)`](@ref rootmeasure). 

`logdensityof` works by first computing `insupport(m, x)`. If this is true, then
`unsafe_logdensityof` is called. If `insupport(m, x)` is known to be `true`, it
can be a little faster to directly call `unsafe_logdensityof(m, x)`. 

To compute log-density relative to `basemeasure(m)` or *define* a log-density
(relative to `basemeasure(m)` or another measure given explicitly), see
`logdensity_def`. 

To compute a log-density relative to a specific base-measure, see
`logdensity_rel`. 
"""
@inline function logdensityof(μ::AbstractMeasure, x)
    result = dynamic(unsafe_logdensityof(μ, x))
    ifelse(insupport(μ, x) == true, result, oftype(result, -Inf))
end

export unsafe_logdensityof

# https://discourse.julialang.org/t/counting-iterations-to-a-type-fixpoint/75876/10?u=cscherrer
"""
    unsafe_logdensityof(m, x)

Compute the log-density of the measure `m` at `x` relative to `rootmeasure(m)`.
This is "unsafe" because it does not check `insupport(m, x)`.

See also `logdensityof`.
"""
@inline function unsafe_logdensityof(μ::M, x) where {M}
    ℓ_0 = logdensity_def(μ, x)
    b_0 = μ
    Base.Cartesian.@nexprs 10 i -> begin  # 10 is just some "big enough" number
        b_{i} = basemeasure(b_{i - 1}, x)
        if b_{i} isa typeof(b_{i - 1})
            return ℓ_{i - 1}
        end
        ℓ_{i} = let Δℓ_{i} = logdensity_def(b_{i}, x)
            ℓ_{i - 1} + Δℓ_{i}
        end
    end
    return ℓ_10
end

export density_rel

@inline density_rel(μ, ν, x) = exp(logdensity_rel(μ, ν, x))

export logdensity_rel

@inline return_type(f, args::Tuple) = Core.Compiler.return_type(f, Tuple{typeof.(args)...})

unstatic(::Type{T}) where {T} = T
unstatic(::Type{StaticFloat64{X}}) where {X} = Float64

"""
    logdensity_rel(m1, m2, x)

Compute the log-density of `m1` relative to `m2` at `x`. This function checks
whether `x` is in the support of `m1` or `m2` (or both, or neither). If `x` is
known to be in the support of both, it can be more efficient to call
`unsafe_logdensity_rel`. 
"""
@inline function logdensity_rel(μ::M, ν::N, x::X) where {M,N,X}
    T = unstatic(
        promote_type(
            return_type(logdensity_def, (μ, x)),
            return_type(logdensity_def, (ν, x)),
        ),
    )
    inμ = insupport(μ, x)
    inν = insupport(ν, x)
    inμ || return convert(T, ifelse(inν, -Inf, NaN))
    inν || return convert(T, Inf)

    return unsafe_logdensity_rel(μ, ν, x)
end

"""
    unsafe_logdensity_rel(m1, m2, x)

Compute the log-density of `m1` relative to `m2` at `x`, assuming `x` is
known to be in the support of both `m1` and `m2`.

See also `logdensity_rel`.
"""
@inline function unsafe_logdensity_rel(μ::M, ν::N, x::X) where {M,N,X}
    if static_hasmethod(logdensity_def, Tuple{M,N,X})
        return logdensity_def(μ, ν, x)
    end
    μs = basemeasure_sequence(μ)
    νs = basemeasure_sequence(ν)
    cb = commonbase(μs, νs, X)
    # _logdensity_rel(μ, ν)
    isnothing(cb) && begin
        μ = μs[end]
        ν = νs[end]
        @warn """
        No common base measure for
            $μ
        and
            $ν

        Returning a relative log-density of NaN. If this is incorrect, add a
        three-argument method
            logdensity_def($μ, $ν, x)
        """
        return NaN
    end
    return _logdensity_rel(μs, νs, cb, x)
end

# Note that this method assumes `μ` and `ν` to have the same type
function logdensity_def(μ::T, ν::T, x) where {T}
    if μ === ν
        return zero(logdensity_def(μ, x))
    else
        return logdensity_def(μ, x) - logdensity_def(ν, x)
    end
end

@generated function _logdensity_rel(
    μs::Tμ,
    νs::Tν,
    ::Tuple{StaticInt{M},StaticInt{N}},
    x::X,
) where {Tμ,Tν,M,N,X}
    sμ = schema(Tμ)
    sν = schema(Tν)

    q = quote
        $(Expr(:meta, :inline))
        ℓ = logdensity_def(μs[$M], νs[$N], x)
    end

    for i in 1:M-1
        push!(q.args, :(Δℓ = logdensity_def(μs[$i], x)))
        # push!(q.args, :(println("Adding", Δℓ)))
        push!(q.args, :(ℓ += Δℓ))
    end

    for j in 1:N-1
        push!(q.args, :(Δℓ = logdensity_def(νs[$j], x)))
        # push!(q.args, :(println("Subtracting", Δℓ)))
        push!(q.args, :(ℓ -= Δℓ))
    end

    push!(q.args, :(return ℓ))
    return q
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
rebase(μ, ν) = ∫(𝒹(μ, ν), ν)
