export PrimitiveMeasure

"""
    abstract type PrimitiveMeasure <: AbstractMeasure end    

In the MeasureTheory ecosystem, a _primitive_ measure is a measure for which the
definition and construction do not depend on any other measure. Primitive
measures satisfy the following laws:

    basemeasure(μ::PrimitiveMeasure) = μ

    logdensity_def(μ::PrimitiveMeasure, x) = 0.0

    logdensity_def(μ::M, ν::M, x) where {M<:PrimitiveMeasure} = 0.0
"""
abstract type PrimitiveMeasure <: AbstractMeasure end

basemeasure(μ::PrimitiveMeasure) = μ

@inline basemeasure_depth(::PrimitiveMeasure) = static(0)

@inline logdensityof(::PrimitiveMeasure, x::Real) = zero(float(typeof(x)))
@inline logdensityof(::PrimitiveMeasure, x) = false

logdensity_def(::PrimitiveMeasure, x) = static(0.0)

logdensity_def(μ::M, ν::M, x) where {M<:PrimitiveMeasure} = 0.0

function Pretty.quoteof(μ::M) where {M<:PrimitiveMeasure}
    :($M())
end
