# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

const DirichletMeasure = AsMeasure{<:Dirichlet}

MeasureBase.getdof(d::Dirichlet) = length(d) - 1
MeasureBase.getdof(m::DirichletMeasure) = getdof(m.obj)
MeasureBase.fast_dof(m::DirichletMeasure) = getdof(m.obj)

@inline MeasureBase.preferred_stdmeasure(::Type{<:Dirichlet}) = StdUniform

# Densities are relative to the coordinate Lebesgue measure on the simplex:
@inline MeasureBase.rootmeasure(m::DirichletMeasure) = Lebesgue(MeasureBase.StdSimplex(length(m.obj)))
