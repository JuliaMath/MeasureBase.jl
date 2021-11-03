export TrivialMeasure

struct TrivialMeasure <: PrimitiveMeasure end

sampletype(::TrivialMeasure) = Nothing

Pretty.quoteof(::TrivialMeasure) = :(TrivialMeasure())
