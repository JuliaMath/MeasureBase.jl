export NormalizedMeasure

"""
    struct NormalizedMeasure{M,L} <: AbstractWeightedMeasure

The measure `base` divided by its mass, see [`normalize`](@ref).

`logmass::L` is the log of the mass of `base`, or an
[`AbstractUnknownMass`](@ref) if it is unknown. Normalized measures of
unknown mass have no density, but support `rand`, transports and
`insupportof`, which only need the normalized law.
"""
struct NormalizedMeasure{M,L} <: AbstractWeightedMeasure
    base::M
    logmass::L
end

"""
    normalize(μ::AbstractMeasure)

The normalized measure `μ / massof(μ)`, a [`NormalizedMeasure`](@ref).

Measures of statically known unit mass (see [`massof`](@ref)), like
standard measures and distributions, are returned as they are. `normalize`
strips weights, `normalize(weightedmeasure(w, μ)) == normalize(μ)`, and is
idempotent. If the mass of `μ` is unknown, the result has no density but
can still be sampled and transported. Throws an `ArgumentError` for
measures of zero or infinite mass.
"""
normalize(μ::AbstractMeasure) = _normalize(μ, massof(μ))
normalize(μ::AbstractWeightedMeasure) = normalize(basemeasure(μ))
normalize(μ::NormalizedMeasure) = μ

@inline _normalize(μ, mass) = NormalizedMeasure(μ, _checked_logmass(μ, log(dynamic(mass))))
@inline _normalize(μ, ::Union{StaticInt{1},StaticFloat64{1.0}}) = μ
@inline _normalize(μ, mass::AbstractUnknownMass) = NormalizedMeasure(μ, mass)

@inline _checked_logmass(μ, ℓ) = (_check_normalizable(μ, isfinite(ℓ)); ℓ)
@inline _check_normalizable(μ, finite::Bool) = finite || _throw_not_normalizable(μ)
# Traced or device values can't be checked:
@inline _check_normalizable(μ, ::Any) = nothing
@noinline _throw_not_normalizable(μ) =
    throw(ArgumentError("Can't normalize a measure of type $(nameof(typeof(μ))), its mass is zero or infinite"))

massof(::NormalizedMeasure) = static(1.0)

_logweight(μ::NormalizedMeasure) = -μ.logmass
_logweight(μ::NormalizedMeasure{<:Any,<:AbstractUnknownMass}) = _throw_unknown_normalization(μ.base)

# Stieltjes measure functions scale with the mass of the base itself
# (not `exp(logmass)`), so that their range is exactly the unit interval:
smf(μ::NormalizedMeasure, x) = _unscaled_smf(_basemass(μ), smf(μ.base, x))
csmf(μ::NormalizedMeasure, x) = _unscaled_smf(_basemass(μ), csmf(μ.base, x))
invsmf(μ::NormalizedMeasure, p) = invsmf(μ.base, p * _basemass(μ))
invcsmf(μ::NormalizedMeasure, p) = invcsmf(μ.base, p * _basemass(μ))

@inline _unscaled_smf(mass, F) = F / mass
@inline _unscaled_smf(mass, F::NoSMF) = F
@inline _basemass(μ::NormalizedMeasure) = massof(μ.base)
@inline _basemass(μ::NormalizedMeasure{<:Any,<:AbstractUnknownMass}) = _throw_unknown_normalization(μ.base)
@noinline _throw_unknown_normalization(μ) =
    throw(ArgumentError("A normalized measure of type $(nameof(typeof(μ))) has no density, its mass is unknown"))

@inline mspace_elsize(μ::NormalizedMeasure) = mspace_elsize(μ.base)
@inline mspace_flatsize(μ::NormalizedMeasure) = mspace_flatsize(μ.base)
@inline mspace_flatsize(::Type{<:NormalizedMeasure{M}}) where {M} = mspace_flatsize(M)
@inline mspace_ndims(::Type{<:NormalizedMeasure{M}}) where {M} = mspace_ndims(M)
@inline stream_ndims(::Type{<:NormalizedMeasure{M}}) where {M} = stream_ndims(M)
@inline fixed_stream_size(::Type{<:NormalizedMeasure{M}}) where {M} = fixed_stream_size(M)
@inline getdof(μ::NormalizedMeasure) = getdof(μ.base)
@inline fast_dof(μ::NormalizedMeasure) = fast_dof(μ.base)
@inline preferred_stdmeasure(::Type{<:NormalizedMeasure{M}}) where {M} = preferred_stdmeasure(M)

insupportof(μ::NormalizedMeasure, x) = insupportof(μ.base, x)
gentype(μ::NormalizedMeasure) = gentype(μ.base)

Pretty.tile(μ::NormalizedMeasure) = Pretty.list_layout([Pretty.tile(μ.base)]; prefix = :normalize)

Adapt.adapt_structure(to, μ::NormalizedMeasure) = NormalizedMeasure(Adapt.adapt(to, μ.base), μ.logmass)

# Transport Jacobians are density ratios of the laws that transports map:
# measures without their weights, structurally (through weighted measures,
# tuple products, powers, combinations, pushforwards and binds), normalized
# at the leaves. Leaves of unknown mass have no normalized density:
@inline _unweighted(μ) = _normalized_leaf(μ, massof(μ))
@inline _normalized_leaf(μ, mass) = _normalize(μ, mass)
@noinline _normalized_leaf(μ, ::AbstractUnknownMass) =
    throw(ArgumentError("Transport Jacobians require normalized densities, the mass of measures of type $(nameof(typeof(μ))) is unknown"))
