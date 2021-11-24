export CountingMeasure

struct CountingMeasure <: PrimitiveMeasure end



# gentype(::CountingMeasure{ℝ}) = Float64
# gentype(::CountingMeasure{ℝ₊}) = Float64
# gentype(::CountingMeasure{𝕀}) = Float64

gentype(::CountingMeasure) = Int


logdensity_def(::CountingMeasure, x) = zero(float(x))

# (::CountingMeaure)(s) = length(Set(s))
