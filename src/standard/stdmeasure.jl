abstract type StdMeasure <: AbstractMeasure end

StdMeasure(::typeof(rand)) = StdUniform()
StdMeasure(::typeof(randexp)) = StdExponential()

getdof(::StdMeasure) = static(1)

getdof(μ::PowerMeasure{<:StdMeasure}) = prod(map(length, d.axes))
