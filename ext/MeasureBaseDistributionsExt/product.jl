# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

@static if isdefined(Distributions, :Product)
    MeasureBase.AbstractMeasure(obj::Distributions.Product) = productmeasure(map(asmeasure, obj.v))

    function AsMeasure{D}(::D) where {D<:Distributions.Product}
        throw(ArgumentError("Don't wrap Distributions.Product into MeasureBase.AsMeasure, use asmeasure to convert instead."))
    end
end

function Distributions.Distribution(
    m::MeasureBase.ProductMeasure{<:AbstractArray{<:AsMeasure{<:Distribution{Univariate}}}},
)
    Distributions.product_distribution(map(x -> x.obj, MeasureBase.marginals(m)))
end

function Base.convert(
    ::Type{Distribution},
    m::MeasureBase.ProductMeasure{<:AbstractArray{<:AsMeasure{<:Distribution{Univariate}}}},
)
    Distributions.Distribution(m)
end

@static if isdefined(Distributions, :ProductDistribution)
    MeasureBase.AbstractMeasure(obj::Distributions.ProductDistribution) = productmeasure(map(asmeasure, obj.dists))

    function AsMeasure{D}(::D) where {D<:Distributions.ProductDistribution}
        throw(ArgumentError("Don't wrap Distributions.ProductDistribution into MeasureBase.AsMeasure, use asmeasure to convert instead."))
    end

    function Distributions.Distribution(
        m::MeasureBase.ProductMeasure{<:AbstractArray{<:AsMeasure{<:Distribution}}},
    )
        Distributions.product_distribution(map(x -> x.obj, MeasureBase.marginals(m)))
    end

    function Distributions.Distribution(
        m::MeasureBase.ProductMeasure{<:Tuple{Vararg{AsMeasure{<:Distribution}}}},
    )
        Distributions.product_distribution(map(x -> x.obj, MeasureBase.marginals(m))...)
    end

    function Base.convert(
        ::Type{Distribution},
        m::MeasureBase.ProductMeasure{<:AbstractArray{<:AsMeasure{<:Distribution}}},
    )
        Distributions.Distribution(m)
    end

    function Base.convert(
        ::Type{Distribution},
        m::MeasureBase.ProductMeasure{<:Tuple{Vararg{AsMeasure{<:Distribution}}}},
    )
        Distributions.Distribution(m)
    end
end
