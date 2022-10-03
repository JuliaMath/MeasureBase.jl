export TrivialMeasure

struct TrivialMeasure <: PrimitiveMeasure end

gentype(::TrivialMeasure) = Nothing

insupport(::TrivialMeasure, x) = False

massof(::TrivialMeasure) = static(0.0)