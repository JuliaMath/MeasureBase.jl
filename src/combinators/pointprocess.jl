
"""
    struct PointProcessMeasure{K,M<:AbstractMeasure} <: AbstractMeasure

A `ν = PointProcessMeasure(Λ)` represents a Poisson random measure, based on an
intensity measure `Λ` that specifies an (non-homogenous) Poisson Point Process
(PPP).

`logdensity_of(ν, x)` requires x to be a vector or set of variates of `ν`.
Vectors are treated as unordered, so a vector `x` is implicitly treated as
`Set(x)`, and both result in the same log-density value.

`PointProcessMeasure(f_kernel, λ::Real)` and
`PointProcessMeasure(f_kernel, λ::Real * Dirac(x)` generate an
[`SingletonPPMeasure`](@ref).
"""
struct PointProcessMeasure{K,M<:AbstractMeasure} <: AbstractMeasure
    f_kernel::K
    Λ::M
    PointProcessMeasure{K,M}(f_kernel::F, Λ::M) where {K,M} = new{K,M}(f_kernel, Λ)
end


pointprocess_measure(f_kernel::F, Λ::M) where {K,M<:AbstractMeasure} = PointProcessMeasure{K,M}(f_kernel, Λ)



function _invalid_intensity_base()
    throw(ArgumentError("When decomposing intensity measure Λ of $(nameof(typeof(ν))) as Λ = λ * μ, the mass of μ must be one."))
end

function _rate_and_dist(ν::PointProcessMeasure{<:WeightedMeasure})
    λ = exp(ν.Λ.logweight)
    μ = ν.Λ.base
    massof(μ) ≈ 1 || _invalid_intensity_base
    return λ, μ
end

function _rate_and_dist(ν::PointProcessMeasure)
    throw(ArgumentError("Don't know how to decompose intensity measure Λ of $(nameof(typeof(ν))) as Λ = λ * μ with a rate λ and a point distribution μ."))
end


insupport(ν::PointProcessMeasure, x::_SomeVectorOrSet) = all(Base.Fix1(insupport, ν.Λ), x)
insupport(ν::PointProcessMeasure, ::Any) = false


function logdensity_def(ν::PointProcessMeasure, x)
    f_kernel = ν.f_kernel
    λ, μ = _rate_and_dist(ν)
    n = length(x)
    Tμ = typeof(μ)
    Tx_i = eltype(x)
    R = Core.Compiler.return_type(logdensity_def, Tuple{Tμ,Tx_i})
    ld_totalrate = R(logdensity_def(f_kernel(λ), n))::R
    ld_points = isempty(x) ? zero(R) : R(sum(Base.Fix1(logdensity_def, μ), x))::R
    return ld_totalrate + ld_points
end

@inline basemeasure(::PointProcessMeasure) = PowersUnionMeasure(basemeasure(ν.Λ))
@inline rootmeasure(::PointProcessMeasure) = PowersUnionMeasure(rootmeasure(ν.Λ))


@inline getdof(ν::PointProcessMeasure) = NoDOF(typeof(ν))

function Base.rand(rng::AbstractRNG, ::Type{T}, ν::PointProcessMeasure) where {T}
    f_kernel = ν.f_kernel
    λ, μ = _rate_and_dist(ν)
    n = rand(rng, f_kernel(λ))
    return rand(rng, T, μ^n)
end



"""
    const SingletonPPMeasure{T<:Real} = PointProcessMeasure{F,WeightedMeasure{T,Dirac{T}}}

Alias for a [`PointProcessMeasure`](@ref) over a singleton set.

User code should not instantiate `SingletonPPMeasure` directly, use
`pointprocess_measure(f_kernel, λ::Real * Dirac(x))` instead.
"""
const SingletonPPMeasure{F,T<:Real,T} = PointProcessMeasure{F,WeightedMeasure{T,Dirac{T}}}
SingletonPPMeasure(f_kernel, λ::Real) = PointProcessMeasure(f_kernel, λ)

@inline _rate_and_dist(ν::SingletonPPMeasure) = exp(ν.Λ.logweight), ν.Λ.base

@inline insupport(ν::SingletonPPMeasure{T}, ::_SomeVectorOrSet{T}) where T = all(isequal ...) #!!!!!!!!!!

logdensity_def(ν::PointProcessMeasure, x::_SomeVectorOrSet{T}) = logdensity_def(ν, length(x))
@inline insupport(ν::SingletonPPMeasure, ::_SomeVectorOrSet{T}) = true

@inline basemeasure(::SingletonPPMeasure) = Counting()
@inline rootmeasure(::SingletonPPMeasure) = Counting()

function Base.rand(rng::AbstractRNG, ::Type, ν::SingletonPPMeasure)
    f_kernel = ν.f_kernel
    λ, _ = _rate_and_dist(ν)
    n = rand(rng, f_kernel(λ))
    Fill((), n)
end


# # To be added later in the Distributions extension:
# function bind(μ::AsMeasure{Distributions.Poisson}, f_transition::Base.Fix1{typeof(^)})
#     λ, _ = _rate_and_dist(ν)
#     μ = f_transition.x
#     return PointProcessMeasure(Poisson, λ * μ)
# end
