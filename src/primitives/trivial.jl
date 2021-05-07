export TrivialMeasure

struct TrivialMeasure <: AbstractMeasure end

isprimtype(::TrivialMeasure) = true

sampletype(::TrivialMeasure) = Nothing
