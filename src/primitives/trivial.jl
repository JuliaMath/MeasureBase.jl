export TrivialMeasure

struct TrivialMeasure <: PrimitiveMeasure end

gentype(::TrivialMeasure) = Nothing

insupport(::TrivialMeasure, x) = static(false)

logdensityof_impl(::TrivialMeasure, x) = _neg_inf_logd(x)

massof(::TrivialMeasure) = static(0.0)
