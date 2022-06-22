function superpose(::T, ::T) where {T<:SuperpositionMeasure}
    @error "FIXME"
end

function kernel(::Type{M}, ::NamedTuple{()}) where {M<:ParameterizedMeasure}
    @error "FIXME"
end

function logdensity_def(
    ::T,
    ::S,
    ::Any,
) where {
    T<:(MeasureBase.SuperpositionMeasure{Tuple{A,B}} where {A,B}),
    S<:(MeasureBase.SuperpositionMeasure{Tuple{A,B}} where {A,B}),
}
    @error "FIXME"
end

function transport_def(::StdUniform, ::StdExponential, ::NoTransformOrigin)
    @error "FIXME"
end

@inline function transport_def(::StdExponential, ::StdUniform, ::NoTransformOrigin)
    @error "FIXME"
end

function transport_def(::StdUniform, ::StdExponential, ::NoTransport)
    @error "FIXME"
end

@inline function transport_def(::StdExponential, ::StdUniform, ::NoTransport)
    @error "FIXME"
end

function transport_def(::StdUniform, ::StdLogistic, ::NoTransformOrigin)
    @error "FIXME"
end

function transport_def(::StdUniform, ::StdLogistic, ::NoTransport)
    @error "FIXME"
end

function transport_def(::StdLogistic, ::StdUniform, ::NoTransport)
    @error "FIXME"
end

@inline function transport_def(::StdLogistic, ::StdUniform, ::NoTransformOrigin)
    @error "FIXME"
end

@inline function transport_def(::MU, ::MU, ::NoTransport) where {MU<:StdMeasure}
    @error "FIXME"
end

@inline function transport_def(::MU, ::MU, ::NoTransformOrigin) where {MU<:StdMeasure}
    @error "FIXME"
end

function transport_def(ν::StdMeasure, μ::PowerMeasure{<:StdMeasure}, x)
    return transport_def(ν, μ.parent, only(x))
end

function transport_def(::StdMeasure, ::PowerMeasure{<:StdMeasure}, ::NoTransport)
    @error "FIXME"
end

function transport_def(::StdMeasure, ::PowerMeasure{<:StdMeasure}, ::NoTransformOrigin)
    @error "FIXME"
end

function transport_def(ν::PowerMeasure{<:StdMeasure}, μ::StdMeasure, x)
    return Fill(transport_def(ν.parent, μ, only(x)), map(length, ν.axes)...)
end

function transport_def(::PowerMeasure{<:StdMeasure}, ::StdMeasure, ::NoTransformOrigin)
    @error "FIXME"
end

function transport_def(::PowerMeasure{<:StdMeasure}, ::StdMeasure, ::NoTransport)
    @error "FIXME"
end

function transport_def(
    ::PowerMeasure{<:StdMeasure,<:Tuple{Base.OneTo}},
    ::PowerMeasure{<:StdMeasure,<:Tuple{Base.OneTo}},
    ::NoTransport,
)
    @error "FIXME"
end

function transport_def(
    ::PowerMeasure{<:StdMeasure,<:Tuple{Base.OneTo}},
    ::PowerMeasure{<:StdMeasure,<:Tuple{Base.OneTo}},
    ::NoTransformOrigin,
)
    @error "FIXME"
end

function transport_def(
    ::PowerMeasure{<:StdMeasure,<:Tuple{Vararg{Base.OneTo,N}}},
    ::PowerMeasure{<:StdMeasure,<:Tuple{Vararg{Base.OneTo,M}}},
    ::NoTransport,
) where {N,M}
    @error "FIXME"
end

function transport_def(
    ::PowerMeasure{<:StdMeasure,<:Tuple{Vararg{Base.OneTo,N}}},
    ::PowerMeasure{<:StdMeasure,<:Tuple{Vararg{Base.OneTo,M}}},
    ::NoTransformOrigin,
) where {N,M}
    @error "FIXME"
end

function transport_to(::Type{NU}, ::Type{MU}) where {MU<:StdMeasure,NU<:StdMeasure}
    @error "FIXME"
end

function transport_def(::Dirac, ::PowerMeasure{<:StdMeasure}, ::NoTransport)
    @error "FIXME"
end

function transport_def(::Dirac, ::PowerMeasure{<:StdMeasure}, ::NoTransformOrigin)
    @error "FIXME"
end

@inline function transport_def(::PowerMeasure{<:StdMeasure}, ::Dirac, ::NoTransport)
    @error "FIXME"
end

@inline function transport_def(::PowerMeasure{<:StdMeasure}, ::Dirac, ::NoTransformOrigin)
    @error "FIXME"
end
