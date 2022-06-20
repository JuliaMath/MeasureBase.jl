abstract type StdMeasure <: AbstractMeasure end

StdMeasure(::typeof(rand)) = StdUniform()
StdMeasure(::typeof(randexp)) = StdExponential()


@inline check_dof(::StdMeasure, ::StdMeasure) = nothing


@inline transport_def(::MU, μ::MU, x) where {MU<:StdMeasure} = x

@inline function transport_def(::MU, ::MU, ::NoTransport) where MU<:StdMeasure
    @error "FIXME"
end

@inline function transport_def(::MU, ::MU, ::NoTransformOrigin) where MU<:StdMeasure
    @error "FIXME"
end

function transport_def(ν::StdMeasure, μ::PowerMeasure{<:StdMeasure}, x)
    return transport_def(ν, μ.parent, only(x))
end

transport_def(::MeasureBase.StdMeasure, ::MeasureBase.PowerMeasure{<:MeasureBase.StdMeasure}, ::MeasureBase.NoTransport) = @error "FIXME"

transport_def(::MeasureBase.StdMeasure, ::MeasureBase.PowerMeasure{<:MeasureBase.StdMeasure}, ::MeasureBase.NoTransformOrigin) = @error "FIXME"

function transport_def(ν::PowerMeasure{<:StdMeasure}, μ::StdMeasure, x)
    return Fill(transport_def(ν.parent, μ, only(x)), map(length, ν.axes)...)
end

transport_def(::MeasureBase.PowerMeasure{<:MeasureBase.StdMeasure}, ::MeasureBase.StdMeasure, ::MeasureBase.NoTransformOrigin) = @error "FIXME"

function transport_def(::PowerMeasure{<:StdMeasure}, ::StdMeasure, ::NoTransport)
    @error "FIXME"
end

function transport_def(ν::PowerMeasure{<:StdMeasure,<:NTuple{1,Base.OneTo}}, μ::PowerMeasure{<:StdMeasure,<:NTuple{1,Base.OneTo}}, x)
    return transport_to(ν.parent, μ.parent).(x)
end

transport_def(::MeasureBase.PowerMeasure{<:MeasureBase.StdMeasure, <:Tuple{Base.OneTo}}, ::MeasureBase.PowerMeasure{<:MeasureBase.StdMeasure, <:Tuple{Base.OneTo}}, ::MeasureBase.NoTransport) = @error "FIXME"

function transport_def(::PowerMeasure{<:StdMeasure, <:Tuple{Base.OneTo}}, ::PowerMeasure{<:StdMeasure, <:Tuple{Base.OneTo}}, ::NoTransformOrigin)
    @error "FIXME"
end

function transport_def(ν::PowerMeasure{<:StdMeasure,<:NTuple{N,Base.OneTo}}, μ::PowerMeasure{<:StdMeasure,<:NTuple{M,Base.OneTo}}, x) where {N,M}
    return reshape(transport_to(ν.parent, μ.parent).(x), map(length, ν.axes)...)
end

function transport_def(::MeasureBase.PowerMeasure{<:MeasureBase.StdMeasure, <:Tuple{Vararg{Base.OneTo, N}}}, ::MeasureBase.PowerMeasure{<:MeasureBase.StdMeasure, <:Tuple{Vararg{Base.OneTo, M}}}, ::MeasureBase.NoTransport) where {N, M}
    @error "FIXME"
end

function transport_def(::PowerMeasure{<:StdMeasure, <:Tuple{Vararg{Base.OneTo, N}}}, ::PowerMeasure{<:StdMeasure, <:Tuple{Vararg{Base.OneTo, M}}}, ::NoTransformOrigin) where {N, M}
    @error "FIXME"
end

# Implement transport_to(NU::Type{<:StdMeasure}, μ) and transport_to(ν, MU::Type{<:StdMeasure}):

_std_measure(::Type{M}, ::StaticInt{1}) where {M<:StdMeasure} = M()
_std_measure(::Type{M}, dof::Integer) where {M<:StdMeasure} = M()^dof
_std_measure_for(::Type{M}, μ::Any) where {M<:StdMeasure} = _std_measure(M, getdof(μ))

transport_to(::Type{NU}, μ) where {NU<:StdMeasure} = transport_to(_std_measure_for(NU, μ), μ)
transport_to(ν, ::Type{MU}) where {MU<:StdMeasure} = transport_to(ν, _std_measure_for(MU, ν))

transport_to(::Type{NU}, ::Type{MU}) where {MU<:MeasureBase.StdMeasure, NU<:MeasureBase.StdMeasure} = @error "FIXME"

# Transform between standard measures and Dirac:

@inline transport_def(ν::Dirac, ::PowerMeasure{<:StdMeasure}, ::Any) = ν.x

transport_def(::MeasureBase.Dirac, ::MeasureBase.PowerMeasure{<:MeasureBase.StdMeasure}, ::MeasureBase.NoTransport) = @error "FIXME"

transport_def(::MeasureBase.Dirac, ::MeasureBase.PowerMeasure{<:MeasureBase.StdMeasure}, ::MeasureBase.NoTransformOrigin) = @error "FIXME"

@inline function transport_def(ν::PowerMeasure{<:StdMeasure}, ::Dirac, ::Any)
    Zeros{Bool}(map(_ -> 0, ν.axes))
end

@inline function transport_def(::PowerMeasure{<:StdMeasure}, ::Dirac, ::NoTransport)
    @error "FIXME"
end

@inline function transport_def(::PowerMeasure{<:StdMeasure}, ::Dirac, ::NoTransformOrigin)
    @error "FIXME"
end