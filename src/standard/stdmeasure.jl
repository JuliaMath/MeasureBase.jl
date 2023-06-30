abstract type StdMeasure <: AbstractMeasure end


const _PowerStdMeasure{N,MU<:StdMeasure} = PowerMeasure{MU,<:NTuple{N,Base.OneTo}}

_get_inner_stdmeasure(::_PowerStdMeasure{N,MU}) where {N,MU} = M()


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

function transport_def(ν::_PowerStdMeasure{MU,1}, μ::_PowerStdMeasure{NU,1}, x,) where {MU,NU}
    return transport_to(ν.parent, μ.parent).(x)
end

transport_origin(μ::_PowerStdMeasure{N}) = ν.parent^product(map(length, μ.axes))

function from_origin(μ::_PowerStdMeasure{N}, x_origin::AbstractVector{<:Real})
    return reshape(x_origin, map(length, μ.axes)...)
end


# Transport to a multivariate standard measure from any measure:

function transport_def(ν::_PowerStdMeasure{1}, μ::AbstractMeasure, ab)
    ν_inner = _get_inner_stdmeasure(ν)
    transport_to_mvstd(ν_inner, μ, ab)
end

function transport_to_mvstd(ν_inner::StdMeasure, μ::AbstractMeasure, x)
    return _to_mvstd_withdof(ν_inner, μ, getdof(μ), x, origin)
end

function _to_mvstd_withdof(ν_inner::StdMeasure, μ::AbstractMeasure, dof_μ, x)
    y = transport_to(ν_inner^dof_μ, μ, x)
    return y
end

function _to_mvstd_withdof(ν_inner::StdMeasure, μ::AbstractMeasure, ::NoDOF, x)
    _to_mvstd_withorigin(ν_inner, μ, transport_origin(μ), x)
end

function _to_mvstd_withorigin(ν_inner::StdMeasure, ::AbstractMeasure, μ_origin, x)
    x_origin = transport_to_mvstd(ν_inner, μ_origin, x)
    from_origin(x_origin)
end

function _to_mvstd_withorigin(ν_inner::StdMeasure, μ::AbstractMeasure, NoTransportOrigin, x)
    throw(ArgumentError("Don't know how to transport values of type $(nameof(typeof(x))) from $(nameof(typeof(μ))) to a power of $(nameof(typeof(ν_inner)))"))
end


# Transport from a multivariate standard measure to any measure:

function transport_def(ν::AbstractMeasure, μ::_PowerStdMeasure{1}, x)
    μ_inner = _get_inner_stdmeasure(μ)
    _transport_from_mvstd(ν, μ_inner, x)
end

function _transport_from_mvstd(ν::AbstractMeasure, μ_inner::StdMeasure, x)
    # Sanity check, should be checked by transport machinery already:
    @assert getdof(μ) == length(eachindex(x)) && x isa AbstractVector
    y, x_rest = transport_from_mvstd_with_rest(ν, μ_inner, x)
    if !isempty(x_rest)
        throw(ArgumentError("Input value too long during transport"))
    end
    return y
end

function transport_from_mvstd_with_rest(ν::AbstractMeasure, μ_inner::StdMeasure, x)
    dof_ν = getdof(ν)
    origin = transport_origin(ν)
    return _from_mvstd_with_rest_withdof(ν, getdof(ν), μ_inner, x, dof_ν, origin)
end

function _from_mvstd_with_rest_withdof(ν::AbstractMeasure, dof_ν, μ_inner::StdMeasure, x)
    len_x = length(eachindex(x))

    # Since we can't check DOF of original Bind, we could "run out x" if
    # the original x was too short. `transport_to` below will detect this, but better
    # throw a more informative exception here:
    if len_x < dof_ν
        throw(ArgumentError("Variate too short during transport involving Bind"))
    end

    x_inner_dof, x_rest = _split_after(x, dof_ν)
    y = transport_to(ν, μ_inner^dof_ν, x_inner_dof)
    return y, x_rest
end

function _from_mvstd_with_rest_withdof(ν::AbstractMeasure, ::NoDOF, μ_inner::StdMeasure, x)
    _from_mvstd_with_rest_withorigin(ν, transport_origin(ν), μ_inner, x)
end

function _from_mvstd_with_rest_withorigin(::AbstractMeasure, ν_origin, μ_inner::StdMeasure, x)
    x_origin, x_rest = transport_from_mvstd_with_rest(ν_origin, x, μ_inner)
    from_origin(x_origin), x_rest
end

function _from_mvstd_with_rest_withorigin(ν::AbstractMeasure, NoTransportOrigin, μ_inner::StdMeasure, x)
    throw(ArgumentError("Don't know how to transport value of type $(nameof(typeof(x))) from power of $(nameof(typeof(μ_inner))) to $(nameof(typeof(ν)))"))
end


_empty_zero(::AbstractVector{T}) where {T<:Real} = Fill(zero(T), 0)


# Implement transport_to(NU::Type{<:StdMeasure}, μ) and transport_to(ν, MU::Type{<:StdMeasure})
# for user convenience:

# ToDo: Handle combined/bind measures that don't have a fast getdof!

_std_measure(::Type{M}, ::StaticInteger{1}) where {M<:StdMeasure} = M()
_std_measure(::Type{M}, dof::IntegerLike) where {M<:StdMeasure} = M()^dof
_std_measure_for(::Type{M}, μ::Any) where {M<:StdMeasure} = _std_measure(M, getdof(μ))

function transport_to(ν, ::Type{MU}) where {MU<:StdMeasure}
    transport_to(ν, _std_measure_for(MU, ν))
end

function transport_to(::Type{NU}, μ) where {NU<:StdMeasure}
    transport_to(_std_measure_for(NU, μ), μ)
end



# Transform between standard measures and Dirac:

@inline transport_from_mvstd_with_rest(ν::Dirac, ::StdMeasure, x::Any) = ν.x, x

@inline transport_to_mvstd(ν::PowerMeasure{<:StdMeasure}, ::Dirac, ::Any) = Zeros{Bool}(map(_ -> 0, ν.axes))



#!!!!!!!!!!!!!!!!!!!!!! TODO:

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
