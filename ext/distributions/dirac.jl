# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

MeasureBase.AbstractMeasure(obj::Distributions.Dirac) = MeasureBase.Dirac(obj.value)

function AsMeasure{D}(::D) where {D<:Distributions.Dirac}
    throw(ArgumentError("Don't wrap Distributions.Dirac into MeasureBase.AsMeasure, use asmeasure to convert instead."))
end


Distributions.Distribution(m::MeasureBase.Dirac{<:Real}) = Distribtions.Dirac(m.x)

function Distributions.Distribution(@nospecialize(m::MeasureBase.Dirac{T})) where T
    throw(ArgumentError("Can only convert MeasureBase.Dirac{<:Real} to Distributions.Dirac, but not MeasureBase.Dirac{<:$(nameof(T))}"))
end
