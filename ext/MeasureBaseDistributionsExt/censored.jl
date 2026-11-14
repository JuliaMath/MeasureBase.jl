# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Censored continuous distributions have atoms at their finite bounds. They
# become superpositions of the truncation of the uncensored distribution to
# the bounds and of point masses at the bounds that carry the tail masses
# beyond them. Missing and infinite bounds have no atoms:
const _ContinuousCensored = Distributions.Censored{<:Distribution{Univariate,Continuous},Continuous}

function MeasureBase.AbstractMeasure(d::_ContinuousCensored)
    d0 = d.uncensored
    lo, hi = _censoring_bound(d.lower), _censoring_bound(d.upper)
    t = truncate(asmeasure(d0), MeasureBase.IntervalSets.ClosedInterval(something(lo, -Inf), something(hi, Inf)))
    superpose(t, _bound_atoms(MeasureBase.logsmf, d0, lo)..., _bound_atoms(MeasureBase.logcsmf, d0, hi)...)
end

@inline _censoring_bound(x::Real) = isfinite(x) ? x : nothing
@inline _censoring_bound(::Nothing) = nothing

@inline _bound_atoms(::Any, ::Distribution, ::Nothing) = ()
@inline _bound_atoms(logtail::F, d0::Distribution, x::Real) where {F} = (MeasureBase.weightedmeasure(logtail(d0, x), MeasureBase.Dirac(x)),)

const _CensoredMeasure = SuperpositionMeasure{
    <:Tuple{MeasureBase.TruncatedMeasure{<:DistributionMeasure{Univariate,Continuous}},Vararg{WeightedMeasure{<:Any,<:MeasureBase.Dirac}}},
}

# Back to a censored distribution where the atoms lie at the bounds of the
# truncation and carry the tail masses of the uncensored distribution:
function Distributions.Distribution(m::_CensoredMeasure)
    t, atoms = first(m.components), Base.tail(m.components)
    d0 = Distribution(t.parent)
    at_lo = _bound_atom(atoms, t.lo, MeasureBase.logsmf(d0, t.lo))
    at_hi = _bound_atom(atoms, t.hi, MeasureBase.logcsmf(d0, t.hi))
    count(!isnothing, (at_lo, at_hi)) == length(atoms) || throw(ArgumentError(
        "Only superpositions of a truncation and of point masses with the tail masses at its bounds convert to censored distributions"
    ))
    Distributions.censored(d0, at_lo, at_hi)
end

function _bound_atom(atoms::Tuple, x::Real, logmass::Real)
    any(a -> a.base.x == x && isapprox(exp(a.logweight), exp(logmass)), atoms) ? x : nothing
end

Base.convert(::Type{Distribution}, m::_CensoredMeasure) = Distributions.Distribution(m)
