# Direct transports between standard measures, tail-accurate in both
# directions. Transports via StdUniform lose the upper tail of unbounded
# measures, since the uniform variate saturates at one.

# Standard normal log-cdf and log-ccdf:
@inline _normlogcdf(z) = logerfc(-z * invsqrt2) - logtwo
@inline _normlogccdf(z) = logerfc(z * invsqrt2) - logtwo

# Complementary standard normal cdf, accurate for large positive arguments:
@inline _normccdf(z) = erfc(z * invsqrt2) / 2

@inline function transport_def(::StdExponential, ::StdNormal, z)
    ifelse(z < zero(z), -log1p(-Φ(z)), -log(_normccdf(z)))
end

@inline function transport_def(::StdNormal, ::StdExponential, x)
    ifelse(x < oftype(x, logtwo), Φinv(-expm1(-x)), -Φinv(exp(-x)))
end

@inline transport_def(::StdLogistic, ::StdNormal, z) = _normlogcdf(z) - _normlogccdf(z)

@inline function transport_def(::StdNormal, ::StdLogistic, l)
    ifelse(l < zero(l), Φinv(logistic(l)), -Φinv(logistic(-l)))
end

@inline transport_def(::StdLogistic, ::StdExponential, x) = log(-expm1(-x)) + x

@inline transport_def(::StdExponential, ::StdLogistic, l) = log1pexp(l)


"""
    MeasureBase.stdconvert(::Type{S}, ::Type{T}, x)

Convert a variate `x` of the standard measure type `T` into a variate of
the standard measure type `S`, elementwise for arrays.
"""
function stdconvert end

@inline stdconvert(::Type{S}, ::Type{S}, x) where {S<:StdMeasure} = x
@inline stdconvert(::Type{S}, ::Type{T}, x) where {S<:StdMeasure,T<:StdMeasure} = _StdConvert{S,T}()(x)

struct _StdConvert{S,T} <: Function end
@inline (::_StdConvert{S,T})(x::Number) where {S,T} = transport_def(S(), T(), x)
@inline (k::_StdConvert)(x::AbstractArray) = broadcast(k, x)


"""
    MeasureBase.StdPowerMeasure{MU<:StdMeasure,N}

The type of an `N`-dimensional power of a standard measure of type `MU`.
"""
const StdPowerMeasure{MU<:StdMeasure,N} = PowerMeasure{MU,<:NTuple{N,OneToLike}}

# Powers of standard measures transport directly, by elementwise conversion:
function transport_def(ν::StdPowerMeasure{NU}, μ::StdPowerMeasure{MU}, x) where {NU<:StdMeasure,MU<:StdMeasure}
    _pwr_variate(ν, maybestatic_reshape(stdconvert(NU, MU, x), mspace_flatsize(ν)))
end
