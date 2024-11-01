# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

function MeasureBase.AbstractMeasure(d::Distributions.ReshapedDistribution)
    orig_dist = d.dist
    pushfwd(Reshape(size(d), size(orig_dist)), AbstractMeasure(orig_dist))
end

function AsMeasure{D}(::D) where {D<:Distributions.ReshapedDistribution}
    throw(ArgumentError("Don't wrap Distributions.ReshapedDistribution into MeasureBase.AsMeasure, use asmeasure to convert instead."))
end


# ToDo: Conversion back to Distribution
