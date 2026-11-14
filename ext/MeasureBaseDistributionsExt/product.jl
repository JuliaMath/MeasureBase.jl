# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

@static if isdefined(Distributions, :Product)
    MeasureBase.AbstractMeasure(obj::Distributions.Product) = productmeasure(map(asmeasure, obj.v))

    @inline function MeasureBase.preferred_stdmeasure(::Type{<:Distributions.Product{<:Any,T}}) where {T}
        MeasureBase.preferred_stdmeasure(T)
    end

    function AsMeasure{D}(::D) where {D<:Distributions.Product}
        throw(ArgumentError("Don't wrap Distributions.Product into MeasureBase.AsMeasure, use asmeasure to convert instead."))
    end
end

# Products over arrays of distributions convert to product measures over
# struct arrays, which isn't type-stable. The bare products keep type-stable
# forms of the base measures and test values of these:
_marginal_dists(d::Distributions.Product) = d.v
@static if isdefined(Distributions, :ProductDistribution)
    _marginal_dists(d::Distributions.ProductDistribution{N,M,<:AbstractArray}) where {N,M} = d.dists
    const _ArrayProductDistribution = Union{Distributions.Product,Distributions.ProductDistribution{<:Any,<:Any,<:AbstractArray}}
else
    const _ArrayProductDistribution = Distributions.Product
end

MeasureBase.basemeasure(d::_ArrayProductDistribution) = MeasureBase.rootmeasure(d)
MeasureBase.rootmeasure(d::_ArrayProductDistribution) = _product_root(d, Distributions.value_support(typeof(d)))
_product_root(d, ::Type{Continuous}) = Lebesgue()^size(d)
_product_root(d, ::Type{Discrete}) = Counting()^size(d)
MeasureBase.testvalue(d::_ArrayProductDistribution) = map(testvalue, _marginal_dists(d))
MeasureBase.testvalue(::Type{T}, d::_ArrayProductDistribution) where {T} = map(c -> testvalue(T, c), _marginal_dists(d))

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

    @inline function MeasureBase.preferred_stdmeasure(::Type{<:Distributions.ProductDistribution{N,M,D}}) where {N,M,D<:AbstractArray}
        MeasureBase.preferred_stdmeasure(eltype(D))
    end

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
