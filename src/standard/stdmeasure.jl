abstract type StdMeasure <: AbstractMeasure end


const _PowerStdMeasure{N,M<:StdMeasure} = PowerMeasure{M,<:NTuple{N,Base.OneTo}}

_get_inner_stdmeasure(μ::_PowerStdMeasure{N,M}) where {N,M} = M()


StdMeasure(::typeof(rand)) = StdUniform()
StdMeasure(::typeof(randexp)) = StdExponential()
StdMeasure(::typeof(randn)) = StdNormal()

@inline check_dof(::StdMeasure, ::StdMeasure) = nothing

@inline transport_def(::MU, μ::MU, x) where {MU<:StdMeasure} = x

function transport_def(ν::StdMeasure, μ::PowerMeasure{<:StdMeasure}, x)
    return transport_def(ν, μ.parent, only(x))
end

function transport_def(ν::PowerMeasure{<:StdMeasure}, μ::StdMeasure, x)
    return fill_with(transport_def(ν.parent, μ, only(x)), map(length, ν.axes))
end

function transport_def(
    ν::PowerMeasure{<:StdMeasure,<:NTuple{1,Base.OneTo}},
    μ::PowerMeasure{<:StdMeasure,<:NTuple{1,Base.OneTo}},
    x,
)
    return transport_to(ν.parent, μ.parent).(x)
end

function transport_def(
    ν::PowerMeasure{<:StdMeasure,<:NTuple{N,Base.OneTo}},
    μ::PowerMeasure{<:StdMeasure,<:NTuple{M,Base.OneTo}},
    x,
) where {N,M}
    return reshape(transport_to(ν.parent, μ.parent).(x), map(length, ν.axes)...)
end

# Implement transport_to(NU::Type{<:StdMeasure}, μ) and transport_to(ν, MU::Type{<:StdMeasure}):

_std_measure(::Type{M}, ::StaticInteger{1}) where {M<:StdMeasure} = M()
_std_measure(::Type{M}, dof::IntegerLike) where {M<:StdMeasure} = M()^dof
_std_measure_for(::Type{M}, μ::Any) where {M<:StdMeasure} = _std_measure(M, getdof(μ))

function transport_to(::Type{NU}, μ) where {NU<:StdMeasure}
    transport_to(_std_measure_for(NU, μ), μ)
end

function transport_to(ν, ::Type{MU}) where {MU<:StdMeasure}
    transport_to(ν, _std_measure_for(MU, ν))
end

# Transform between standard measures and Dirac:

@inline transport_def(ν::Dirac, ::PowerMeasure{<:StdMeasure}, ::Any) = ν.x

@inline function transport_def(ν::PowerMeasure{<:StdMeasure}, ::Dirac, ::Any)
    Zeros{Bool}(map(_ -> 0, ν.axes))
end

# Helpers for product transforms and similar:

struct _TransportToStd{NU<:StdMeasure} <: Function end
_TransportToStd{NU}(μ, x) where {NU} = transport_to(NU()^getdof(μ), μ)(x)

struct _TransportFromStd{MU<:StdMeasure} <: Function end
_TransportFromStd{MU}(ν, x) where {MU} = transport_to(ν, MU()^getdof(ν))(x)

function _tuple_transport_def(
    ν::PowerMeasure{NU},
    μs::Tuple,
    xs::Tuple,
) where {NU<:StdMeasure}
    reshape(vcat(map(_TransportToStd{NU}, μs, xs)...), ν.axes)
end

function transport_def(
    ν::PowerMeasure{NU},
    μ::ProductMeasure{<:Tuple},
    x,
) where {NU<:StdMeasure}
    _tuple_transport_def(ν, marginals(μ), x)
end

function transport_def(
    ν::PowerMeasure{NU},
    μ::ProductMeasure{<:NamedTuple{names}},
    x,
) where {NU<:StdMeasure,names}
    _tuple_transport_def(ν, values(marginals(μ)), values(x))
end

@inline _offset_cumsum(s, x, y, rest...) = (s, _offset_cumsum(s + x, y, rest...)...)
@inline _offset_cumsum(s, x) = (s,)
@inline _offset_cumsum(s) = ()

function _stdvar_viewranges(μs::Tuple, startidx::IntegerLike)
    N = map(getdof, μs)
    offs = _offset_cumsum(startidx, N...)
    map((o, n) -> o:o+n-1, offs, N)
end

function _tuple_transport_def(
    νs::Tuple,
    μ::PowerMeasure{MU},
    x::AbstractArray{<:Real},
) where {MU<:StdMeasure}
    vrs = _stdvar_viewranges(νs, firstindex(x))
    xs = map(r -> view(x, r), vrs)
    map(_TransportFromStd{MU}, νs, xs)
end

function transport_def(
    ν::ProductMeasure{<:Tuple},
    μ::PowerMeasure{MU},
    x,
) where {MU<:StdMeasure}
    _tuple_transport_def(marginals(ν), μ, x)
end

function transport_def(
    ν::ProductMeasure{<:NamedTuple{names}},
    μ::PowerMeasure{MU},
    x,
) where {MU<:StdMeasure,names}
    NamedTuple{names}(_tuple_transport_def(values(marginals(ν)), μ, x))
end
