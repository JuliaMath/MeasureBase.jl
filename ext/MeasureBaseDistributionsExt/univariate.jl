# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).


@inline MeasureBase.getdof(::Distribution{Univariate}) = static(1)

@inline MeasureBase.check_dof(a::Distribution{Univariate}, b::Distribution{Univariate}) = nothing

_dist_params_numtype(d::Distribution) = real_numtype(typeof(Distributions.params(d)))

@inline function _result_numtype(d::Distribution{Univariate}, x::T) where {T<:Number}
    float(promote_type(T, _dist_params_numtype(d)))
end


# Generic transports between univariate continuous distributions and
# StdLogistic: the log-cdf and log-ccdf keep both tails accurate on the way
# to the standard measure, quantile and complementary quantile on the way
# back. The implementation hooks are specialized for dual numbers in the
# ForwardDiff extension.

@inline MeasureBase.preferred_stdmeasure(::Type{<:Distribution{Univariate,Continuous}}) = StdLogistic

@inline _trafo_logcdf(d::Distribution{Univariate,Continuous}, x::Number) =
    _trafo_logcdf_impl(_dist_params_numtype(d), d, x)
@inline _trafo_logccdf(d::Distribution{Univariate,Continuous}, x::Number) =
    _trafo_logccdf_impl(_dist_params_numtype(d), d, x)
@inline _trafo_quantile(d::Distribution{Univariate,Continuous}, p::Number) =
    _trafo_quantile_impl(_dist_params_numtype(d), d, p)
@inline _trafo_cquantile(d::Distribution{Univariate,Continuous}, p::Number) =
    _trafo_cquantile_impl(_dist_params_numtype(d), d, p)

@inline _trafo_logcdf_impl(::Type{<:Real}, d::Distribution{Univariate,Continuous}, x::Number) =
    Distributions.logcdf(d, x)
@inline _trafo_logccdf_impl(::Type{<:Real}, d::Distribution{Univariate,Continuous}, x::Number) =
    Distributions.logccdf(d, x)
@inline _trafo_quantile_impl(::Type{<:Real}, d::Distribution{Univariate,Continuous}, p::Number) =
    _dist_quantile(d, p)
@inline _trafo_cquantile_impl(::Type{<:Real}, d::Distribution{Univariate,Continuous}, p::Number) =
    _dist_cquantile(d, p)

@inline _dist_quantile(d::Distribution{Univariate,Continuous}, p::Number) = Distributions.quantile(d, p)
@inline _dist_cquantile(d::Distribution{Univariate,Continuous}, p::Number) = Distributions.cquantile(d, p)

# The quantile implementation of Beta only supports Float64:
const _Float64Compatible = Union{Integer,AbstractFloat}
@inline function _dist_quantile(d::Beta{<:_Float64Compatible}, p::_Float64Compatible)
    Distributions.quantile(d, convert(promote_type(Float64, typeof(p)), p))
end
@inline function _dist_cquantile(d::Beta{<:_Float64Compatible}, p::_Float64Compatible)
    Distributions.cquantile(d, convert(promote_type(Float64, typeof(p)), p))
end

# Rounding errors can push quantiles of truncated distributions slightly
# outside of their support:
const _Truncated = Distributions.Truncated{<:Distribution{Univariate,Continuous}}
@inline _dist_quantile(d::_Truncated, p::Real) = _clamp_to_support(d, Distributions.quantile(d, p))
@inline _dist_cquantile(d::_Truncated, p::Real) = _clamp_to_support(d, Distributions.cquantile(d, p))

function _clamp_to_support(d::_Truncated, x::T) where {T<:Real}
    min_x = T(minimum(d))
    max_x = T(maximum(d))
    if x < min_x && isapprox(x, min_x, atol = 4 * eps(T))
        min_x
    elseif x > max_x && isapprox(x, max_x, atol = 4 * eps(T))
        max_x
    else
        x
    end
end

@inline function MeasureBase.transport_to_std(::Type{StdLogistic}, d::Distribution{Univariate,Continuous}, x)
    R = _result_numtype(d, x)
    l = _trafo_logcdf(d, x) - _trafo_logccdf(d, x)
    ifelse(Distributions.insupport(d, x), convert(R, l), convert(R, NaN))
end

@inline function MeasureBase.transport_from_std(::Type{StdLogistic}, d::Distribution{Univariate,Continuous}, l)
    R = _result_numtype(d, l)
    # From the side that keeps the tail:
    x = l < zero(l) ? _trafo_quantile(d, logistic(l)) : _trafo_cquantile(d, logistic(-l))
    convert(R, x)
end


# Location-scale families of standard measures transport by their affine map:

@inline function _affine_to_std(d::Distribution{Univariate}, x::Number)
    z = (x - Distributions.location(d)) / Distributions.scale(d)
    convert(_result_numtype(d, x), z)
end

@inline function _std_to_affine(d::Distribution{Univariate}, z::Number)
    x = muladd(z, Distributions.scale(d), Distributions.location(d))
    convert(_result_numtype(d, z), x)
end

for (D, S) in [
    (Uniform, StdUniform),
    (Logistic, StdLogistic),
    (Normal, StdNormal)
]
    @eval begin
        @inline MeasureBase.preferred_stdmeasure(::Type{<:$D}) = $S
        @inline MeasureBase.transport_to_std(::Type{$S}, d::$D, x) = _affine_to_std(d, x)
        @inline MeasureBase.transport_from_std(::Type{$S}, d::$D, z) = _std_to_affine(d, z)
    end
end

@inline MeasureBase.preferred_stdmeasure(::Type{<:Exponential}) = StdExponential
@inline MeasureBase.transport_to_std(::Type{StdExponential}, d::Exponential, x) =
    convert(_result_numtype(d, x), Distributions.scale(d) \ x)
@inline MeasureBase.transport_from_std(::Type{StdExponential}, d::Exponential, z) =
    convert(_result_numtype(d, z), Distributions.scale(d) * z)


# Affine transformed distributions transport via the underlying distribution:

const _AffineDist = Distributions.AffineDistribution

@inline MeasureBase.preferred_stdmeasure(::Type{<:_AffineDist{<:Any,<:Any,D}}) where {D} =
    MeasureBase.preferred_stdmeasure(D)

@inline function MeasureBase.transport_to_std(::Type{S}, d::_AffineDist, x) where {S<:StdMeasure}
    transport_to_std(S, d.ρ, d.σ \ (x - d.μ))
end
@inline function MeasureBase.transport_from_std(::Type{S}, d::_AffineDist, z) where {S<:StdMeasure}
    muladd(d.σ, transport_from_std(S, d.ρ, z), d.μ)
end
# Disambiguation with the generic univariate transports:
@inline function MeasureBase.transport_to_std(::Type{StdLogistic}, d::_AffineDist, x)
    transport_to_std(StdLogistic, d.ρ, d.σ \ (x - d.μ))
end
@inline function MeasureBase.transport_from_std(::Type{StdLogistic}, d::_AffineDist, z)
    muladd(d.σ, transport_from_std(StdLogistic, d.ρ, z), d.μ)
end
