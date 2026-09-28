export TrivialMeasure

struct TrivialMeasure <: PrimitiveMeasure end

gentype(::TrivialMeasure) = Nothing

insupport(::TrivialMeasure, x) = static(false)

massof(::TrivialMeasure) = static(0.0)
