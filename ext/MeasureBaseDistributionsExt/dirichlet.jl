# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

const DirichletMeasure = AsMeasure{<:Dirichlet}

MeasureBase.getdof(d::Dirichlet) = length(d) - 1
MeasureBase.getdof(m::DirichletMeasure) = getdof(m.obj)

@inline MeasureBase.preferred_stdmeasure(::Type{<:Dirichlet}) = StdUniform
