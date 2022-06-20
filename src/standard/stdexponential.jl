struct StdExponential <: StdMeasure end

export StdExponential

insupport(d::StdExponential, x) = x ≥ zero(x)

@inline logdensity_def(::StdExponential, x) = -x
@inline basemeasure(::StdExponential) = Lebesgue()

@inline transport_def(::StdUniform, μ::StdExponential, x) = -expm1(-x)
@inline transport_def(::StdExponential, μ::StdUniform, x) = -log1p(-x)

function transport_def(
    ::MeasureBase.StdUniform,
    ::MeasureBase.StdExponential,
    ::MeasureBase.NoTransformOrigin,
)
    @error "FIXME"
end

@inline function transport_def(::StdExponential, ::StdUniform, ::NoTransformOrigin)
    @error "FIXME"
end

function transport_def(
    ::MeasureBase.StdUniform,
    ::MeasureBase.StdExponential,
    ::MeasureBase.NoTransport,
)
    @error "FIXME"
end

@inline function transport_def(::StdExponential, ::StdUniform, ::NoTransport)
    @error "FIXME"
end

Base.rand(rng::Random.AbstractRNG, ::Type{T}, ::StdExponential) where {T} = randexp(rng, T)
