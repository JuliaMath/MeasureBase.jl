abstract type PrimitiveMeasure <: AbstractMeasure end

isprimitive(::PrimitiveMeasure) = true
isprimitive(μ) = false

basemeasure(μ::PrimitiveMeasure) = μ
