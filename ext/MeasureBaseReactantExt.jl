# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseReactantExt

using Reactant: Reactant, TracedRNumber
import MeasureBase
using MeasureBase: RealValues, IntegerValues

Base.in(::TracedRNumber{<:Real}, ::RealValues) = true
Base.in(::TracedRNumber{<:Integer}, ::IntegerValues) = true

# CHLO provides erf_inv but no erfc_inv, and `erf_inv(2p - 1)` saturates in
# the tails, so traced code uses Wichura's rational approximation of the
# standard normal quantile (algorithm AS 241, PPND16), branch-free and
# accurate to a few 1e-16 relative in both tails:
MeasureBase.Φinv(p::TracedRNumber{T}) where {T<:Real} = _ppnd16(p, float(T))

const _ppnd16_a = (3.387132872796366608, 133.14166789178437745, 1971.5909503065514427, 13731.693765509461125,
    45921.953931549871457, 67265.770927008700853, 33430.575583588128105, 2509.0809287301226727)
const _ppnd16_b = (1.0, 42.313330701600911252, 687.1870074920579083, 5394.1960214247511077,
    21213.794301586595867, 39307.89580009271061, 28729.085735721942674, 5226.495278852545925)
const _ppnd16_c = (1.42343711074968357734, 4.6303378461565452959, 5.7694972214606914055, 3.64784832476320460504,
    1.27045825245236838258, 0.24178072517745061177, 0.0227238449892691845833, 7.7454501427834140764e-4)
const _ppnd16_d = (1.0, 2.05319162663775882187, 1.6763848301838038494, 0.68976733498510000455,
    0.14810397642748007459, 0.0151986665636164571966, 5.475938084995344946e-4, 1.05075007164441684324e-9)
const _ppnd16_e = (6.6579046435011037772, 5.4637849111641143699, 1.7848265399172913358, 0.29656057182850489123,
    0.026532189526576123093, 0.0012426609473880784386, 2.71155556874348757815e-5, 2.01033439929228813265e-7)
const _ppnd16_f = (1.0, 0.59983220655588793769, 0.13692988092273580531, 0.0148753612908506148525,
    7.868691311456132591e-4, 1.8463183175100546818e-5, 1.4215117583164458887e-7, 2.04426310338993978564e-15)

@inline _ratpoly(x, num, den, ::Type{T}) where {T} = evalpoly(x, map(T, num)) / evalpoly(x, map(T, den))

# All three regions are evaluated on arguments clamped into their range,
# so that the unselected ones stay finite:
function _ppnd16(p, ::Type{T}) where {T<:AbstractFloat}
    q = p - T(0.5)
    q_c = clamp(q, T(-0.425), T(0.425))
    x_c = q_c * _ratpoly(T(0.180625) - q_c * q_c, _ppnd16_a, _ppnd16_b, T)
    r = sqrt(-log(min(p, 1 - p)))
    x_m = _ratpoly(clamp(r, T(1.6), T(5)) - T(1.6), _ppnd16_c, _ppnd16_d, T)
    x_t = _ratpoly(max(r, T(5)) - T(5), _ppnd16_e, _ppnd16_f, T)
    x_tails = ifelse(r <= T(5), x_m, x_t)
    ifelse(abs(q) <= T(0.425), x_c, ifelse(q < 0, -x_tails, x_tails))
end

end # module MeasureBaseReactantExt
