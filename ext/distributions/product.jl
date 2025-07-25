# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

@static if isdefined(Distributions, :Product)
    MeasureBase.AbstractMeasure(obj::Distributions.Product) = productmeasure(map(asmeasure, obj.v))

    function AsMeasure{D}(::D) where {D<:Distributions.Product}
        throw(ArgumentError("Don't wrap Distributions.Product into MeasureBase.AsMeasure, use asmeasure to convert instead."))
    end
end

@static if isdefined(Distributions, :ProductDistribution)
    MeasureBase.AbstractMeasure(obj::Distributions.ProductDistribution) = productmeasure(map(asmeasure, obj.dists))

    function AsMeasure{D}(::D) where {D<:Distributions.ProductDistribution}
        throw(ArgumentError("Don't wrap Distributions.ProductDistribution into MeasureBase.AsMeasure, use asmeasure to convert instead."))
    end
end
