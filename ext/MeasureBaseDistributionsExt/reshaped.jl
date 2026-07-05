# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

MeasureBase.getdof(μ::ReshapedDistribution) = MeasureBase.getdof(μ.dist)

MeasureBase.transport_origin(μ::ReshapedDistribution) = μ.dist

MeasureBase.to_origin(ν::ReshapedDistribution, y) = reshape(y, size(ν.dist))

MeasureBase.from_origin(ν::ReshapedDistribution, x) = reshape(x, ν.dims)


function MeasureBase.AbstractMeasure(d::Distributions.ReshapedDistribution)
    orig_dist = d.dist
    pushfwd(Reshape(size(d), size(orig_dist)), AbstractMeasure(orig_dist))
end

function AsMeasure{D}(::D) where {D<:Distributions.ReshapedDistribution}
    throw(ArgumentError("Don't wrap Distributions.ReshapedDistribution into MeasureBase.AsMeasure, use asmeasure to convert instead."))
end


function Distributions.Distribution(m::PushforwardMeasure{<:Reshape})
    reshape(Distributions.Distribution(m.origin), asnonstatic(m.f.output_size)...)
end

Base.convert(::Type{Distribution}, m::PushforwardMeasure{<:Reshape}) =
    Distributions.Distribution(m)
