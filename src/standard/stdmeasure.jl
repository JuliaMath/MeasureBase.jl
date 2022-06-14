abstract type StdMeasure<:AbstractMeasure end

StdMeasure(::typeof(rand)) = StdUniform()
StdMeasure(::typeod(randn)) = StdNormal()
StdMeasure(::typeof(randexp)) = StdExponential()
