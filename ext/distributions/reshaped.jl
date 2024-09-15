# This file is a part of DistributionMeasures.jl, licensed under the MIT License (MIT).

MeasureBase.getdof(μ::ReshapedDistribution) = MeasureBase.getdof(μ.dist)

MeasureBase.transport_origin(μ::ReshapedDistribution) = μ.dist

MeasureBase.to_origin(ν::ReshapedDistribution, y) = reshape(y, size(ν.dist))

MeasureBase.from_origin(ν::ReshapedDistribution, x) = reshape(x, ν.dims)
