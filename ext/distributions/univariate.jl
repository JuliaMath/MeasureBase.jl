# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).


@inline MeasureBase.getdof(::Distribution{Univariate}) = static(1)

@inline MeasureBase.check_dof(a::Distribution{Univariate}, b::Distribution{Univariate}) = nothing


# Use ForwardDiff for univariate transformations:
@inline function ChainRulesCore.rrule(::typeof(transport_def), ν::Distribution{Univariate}, μ::Distribution{Univariate}, x::Any)
    ChainRulesCore.rrule(fwddiff(transport_def), ν, μ, x)
end
@inline function ChainRulesCore.rrule(::typeof(transport_def), ν::MeasureBase.StdMeasure, μ::Distribution{Univariate}, x::Any)
    ChainRulesCore.rrule(fwddiff(transport_def), ν, μ, x)
end
@inline function ChainRulesCore.rrule(::typeof(transport_def), ν::Distribution{Univariate}, μ::MeasureBase.StdMeasure, x::Any)
    ChainRulesCore.rrule(fwddiff(transport_def), ν, μ, x)
end


# Generic transformations to/from StdUniform via cdf/quantile:


_dist_params_numtype(d::Distribution) = promote_type(map(typeof, Distributions.params(d))...)


@inline _trafo_cdf(d::Distribution{Univariate,Continuous}, x::Real) = _trafo_cdf_impl(_dist_params_numtype(d), d, x)

@inline _trafo_cdf_impl(::Type{<:Real}, d::Distribution{Univariate,Continuous}, x::Real) = Distributions.cdf(d, x)

@inline function _trafo_cdf_impl(::Type{<:Union{Integer,AbstractFloat}}, d::Distribution{Univariate,Continuous}, x::ForwardDiff.Dual{TAG}) where TAG
    x_v = ForwardDiff.value(x)
    u = Distributions.cdf(d, x_v)
    dudx = Distributions.pdf(d, x_v)
    ForwardDiff.Dual{TAG}(u, dudx * ForwardDiff.partials(x))
end


@inline _trafo_quantile(d::Distribution{Univariate,Continuous}, u::Real) = _trafo_quantile_impl(_dist_params_numtype(d), d, u)

@inline _trafo_quantile_impl(::Type{<:Real}, d::Distribution{Univariate,Continuous}, u::Real) = _trafo_quantile_impl_generic(d, u)

@inline function _trafo_quantile_impl(::Type{<:Union{Integer,AbstractFloat}}, d::Distribution{Univariate,Continuous}, u::ForwardDiff.Dual{TAG}) where {TAG}
    x = _trafo_quantile_impl_generic(d, ForwardDiff.value(u))
    dxdu = inv(Distributions.pdf(d, x))
    ForwardDiff.Dual{TAG}(x, dxdu * ForwardDiff.partials(u))
end


@inline _trafo_quantile_impl_generic(d::Distribution{Univariate,Continuous}, u::Real) = Distributions.quantile(d, u)

# Workaround for Beta dist, ForwardDiff doesn't work for parameters:
@inline _trafo_quantile_impl_generic(d::Beta{T}, u::Real) where {T<:ForwardDiff.Dual} = convert(float(typeof(u)), NaN)
# Workaround for Beta dist, current quantile implementation only supports Float64:
@inline function _trafo_quantile_impl_generic(d::Beta{T}, u::Union{Integer,AbstractFloat}) where {T<:Union{Integer,AbstractFloat}}
    Distributions.quantile(d, convert(promote_type(Float64, typeof(u)), u))
end

#=
# ToDo:

# Workaround for rounding errors that can result in quantile values outside of support of Truncated:
@inline function _trafo_quantile_impl_generic(d::Truncated{<:Distribution{Univariate,Continuous}}, u::Real)
    x = Distributions.quantile(d, u)
    T = typeof(x)
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

# Workaround for rounding errors that can result in quantile values outside of support of Truncated:
@inline function _trafo_quantile_impl_generic(d::Truncated{<:Distribution{Univariate,Continuous}}, u::Real)
    x = Distributions.quantile(d, u)
    T = typeof(x)
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
=#


@inline function _result_numtype(d::Distribution{Univariate}, x::T) where {T<:Real}
    float(promote_type(T, eltype(Distributions.params(d))))
    # firsttype(first(typeof(x), promote_type(map(eltype, Distributions.params(d))...)))
end


@inline function MeasureBase.transport_def(::StdUniform, μ::Distribution{Univariate,Continuous}, x)
    R = _result_numtype(μ, x)
    if Distributions.insupport(μ, x)
        y = _trafo_cdf(μ, x)
        convert(R, y)
    else
        convert(R, NaN)
    end
end


@inline function MeasureBase.transport_def(ν::Distribution{Univariate,Continuous}, ::StdUniform, x::T) where T
    R = _result_numtype(ν, x)
    TF = float(T)
    if 0 <= x <= 1
        # Avoid x ≈ 0 and x ≈ 1 to avoid infinite variate values for target distributions with infinite support:
        mod_x = ifelse(x == 0, zero(TF) + eps(TF), ifelse(x == 1, one(TF) - eps(TF), convert(TF, x)))
        y = _trafo_quantile(ν, mod_x)
        convert(R, y)
    else
        convert(R, NaN)
    end
end


# Use standard measures as transformation origin for scaled/translated equivalents:

function _origin_to_affine(ν::Distribution{Univariate}, y::T) where {T<:Real}
    trg_offs, trg_scale = Distributions.location(ν), Distributions.scale(ν)
    x = muladd(y, trg_scale, trg_offs)
    convert(_result_numtype(ν, y), x)
end

function _affine_to_origin(μ::Distribution{Univariate}, x::T) where {T<:Real}
    src_offs, src_scale = Distributions.location(μ), Distributions.scale(μ)
    y = (x - src_offs) / src_scale
    convert(_result_numtype(μ, x), y)
end

for (A, B) in [
    (Uniform, StdUniform),
    (Logistic, StdLogistic),
    (Normal, StdNormal)
]
    @eval begin
        @inline MeasureBase.transport_origin(::$A) = $B()
        @inline MeasureBase.to_origin(ν::$A, y) = _affine_to_origin(ν, y)
        @inline MeasureBase.from_origin(ν::$A, x) = _origin_to_affine(ν, x)
    end
end

@inline MeasureBase.transport_origin(::Exponential) = StdExponential()
@inline MeasureBase.to_origin(ν::Exponential, y) = Distributions.scale(ν) \ y
@inline MeasureBase.from_origin(ν::Exponential, x) = Distributions.scale(ν) * x



# Transform between univariate and single-element power measure

function MeasureBase.transport_def(ν::Distribution{Univariate}, μ::PowerMeasure{<:StdMeasure}, x)
    return transport_def(ν, μ.parent, only(x))
end

function MeasureBase.transport_def(ν::PowerMeasure{<:StdMeasure}, μ::Distribution{Univariate}, x)
    return Fill(transport_def(ν.parent, μ, only(x)), map(length, ν.axes)...)
end


# Transform between univariate and single-element standard multivariate

function MeasureBase.transport_def(ν::Distribution{Univariate}, μ::StandardDist{D,1}, x) where D
    return transport_def(ν, StandardDist{D}(), only(x))
end

function MeasureBase.transport_def(ν::StandardDist{D,1}, μ::Distribution{Univariate}, x) where D
    return Fill(transport_def(StandardDist{D}(), μ, only(x)), size(ν)...)
end
