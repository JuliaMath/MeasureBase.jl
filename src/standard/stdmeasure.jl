abstract type StdMeasure <: AbstractMeasure end

StdMeasure(::typeof(rand)) = StdUniform()
StdMeasure(::typeof(randn)) = StdNormal()
StdMeasure(::typeof(randexp)) = StdExponential()
