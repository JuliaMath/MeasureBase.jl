# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Distributions are normalized, their moments are those of Distributions:
Statistics.mean(m::DistributionMeasure) = Statistics.mean(m.obj)
Statistics.var(m::DistributionMeasure) = Statistics.var(m.obj)
Statistics.cov(m::DistributionMeasure) = Statistics.cov(m.obj)
Statistics.cov(m::DistributionMeasure{Univariate}) = MeasureBase._no_moment(m, :covariance)

function _normal_truncated_moments(t::MeasureBase.TruncatedMeasure{<:AsMeasure{<:Normal}})
    d = t.parent.obj
    μ, σ, lo, hi = promote(d.μ, d.σ, float(t.lo), float(t.hi))
    # Widen before standardizing, which otherwise loses short Float32 intervals.
    μw, σw = μ * 1.0, σ * 1.0
    # Infinite limits are constant with respect to the Normal parameters.
    a = isinf(lo) ? lo * 1.0 : (lo * 1.0 - μw) / σw
    b = isinf(hi) ? hi * 1.0 : (hi * 1.0 - μw) / σw
    m, v = MeasureBase._normal_truncated_moments(a, b)
    return oftype(μ, μw + σw * m), oftype(σ, σw^2 * v)
end

function Statistics.mean(t::MeasureBase.TruncatedMeasure{<:AsMeasure{<:Normal}})
    first(_normal_truncated_moments(t))
end
function Statistics.var(t::MeasureBase.TruncatedMeasure{<:AsMeasure{<:Normal}})
    last(_normal_truncated_moments(t))
end
