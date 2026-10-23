# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Distributions are normalized, their moments are those of Distributions:
Statistics.mean(m::DistributionMeasure) = Statistics.mean(m.obj)
Statistics.var(m::DistributionMeasure) = Statistics.var(m.obj)
Statistics.cov(m::DistributionMeasure) = Statistics.cov(m.obj)
Statistics.cov(m::DistributionMeasure{Univariate}) = MeasureBase._no_moment(m, :covariance)
