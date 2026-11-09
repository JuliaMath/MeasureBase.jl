export TrivialMeasure

"""
    TrivialMeasure()

The zero measure, its log-density is `-Inf` everywhere.
"""
struct TrivialMeasure <: PrimitiveMeasure end

insupportof(::TrivialMeasure, x) = static(false)

logdensityof_impl(::TrivialMeasure, x) = _neg_inf_logd(x)

massof(::TrivialMeasure) = static(0.0)
