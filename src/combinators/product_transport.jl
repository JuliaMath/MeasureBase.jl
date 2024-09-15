# For transport, always pull a PowerMeasure back to one-dimensional PowerMeasure first:

transport_origin(μ::PowerMeasure{<:Any,N}) where N = ν.parent^product(pwr_size(μ))

function from_origin(μ::PowerMeasure{<:Any,N}, x_origin) where N
    # Sanity check, should never fail:
    @assert x_origin isa AbstractVector
    return reshape(x_origin, pwr_size(μ)...)
end


# A one-dimensional PowerMeasure has an origin if it's parent has an origin:

transport_origin(μ::PowerMeasure{<:AbstractMeasure,1}) = _origin_pwr(::typeof(μ), transport_origin(μ.parent), μ.axes)
_pwr_origin(::Type{MU}, parent_origin, axes) = parent_origin^axes
_pwr_origin(::Type{MU}, ::NoTransportOrigin, axes) = NoTransportOrigin{MU}

function from_origin(μ::PowerMeasure{<:AbstractMeasure,1}, x_origin)
    # Sanity check, should never fail:
    @assert x_origin isa AbstractVector
    from_origin.(Ref(μ.parent), x_origin)
end


# Transport between univariate standard measures and power measures of size one:

function transport_def(ν::StdMeasure, μ::PowerMeasure{<:StdMeasure}, x)
    return transport_def(ν, μ.parent, only(x))
end

function transport_def(ν::PowerMeasure{<:StdMeasure}, μ::StdMeasure, x)
    return fill_with(transport_def(ν.parent, μ, only(x)), map(length, ν.axes))
end

function transport_def(ν::StdPowerMeasure{MU,1}, μ::StdPowerMeasure{NU,1}, x,) where {MU,NU}
    return transport_to(ν.parent, μ.parent).(x)
end


# Transport to a multivariate standard measure from any measure:

function transport_def(ν::StdPowerMeasure{MU,1}, μ::AbstractMeasure, ab) where MU
    ν_inner = _inner_stdmeasure(ν)
    transport_to_mvstd(ν_inner, μ, ab)
end

function transport_to_mvstd(ν_inner::StdMeasure, μ::AbstractMeasure, x)
    return _to_mvstd_withdof(ν_inner, μ, fast_dof(μ), x, origin)
end

function _to_mvstd_withdof(ν_inner::StdMeasure, μ::AbstractMeasure, dof_μ::IntegerLike, x)
    y = transport_to(ν_inner^dof_μ, μ, x)
    return y
end

function _to_mvstd_withdof(ν_inner::StdMeasure, μ::AbstractMeasure, ::AbstractNoDOF, x)
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

function transport_def(ν::AbstractMeasure, μ::StdPowerMeasure{MU,1}, x) where MU
    μ_inner = _inner_stdmeasure(μ)
    _transport_from_mvstd(ν, μ_inner, x)
end

function _transport_from_mvstd(ν::AbstractMeasure, μ_inner::StdMeasure, x)
    y, x_rest = transport_from_mvstd_with_rest(ν, μ_inner, x)
    if !isempty(x_rest)
        throw(ArgumentError("Input value too long during transport"))
    end
    return y
end

function transport_from_mvstd_with_rest(ν::AbstractMeasure, μ_inner::StdMeasure, x)
    dof_ν = fast_dof(ν)
    origin = transport_origin(ν)
    return _from_mvstd_with_rest_withdof(ν, dof_ν, μ_inner, x, dof_ν, origin)
end

function _from_mvstd_with_rest_withdof(ν::AbstractMeasure, dof_ν::IntegerLike, μ_inner::StdMeasure, x)
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

function _from_mvstd_with_rest_withdof(ν::AbstractMeasure, ::AbstractNoDOF, μ_inner::StdMeasure, x)
    _from_mvstd_with_rest_withorigin(ν, transport_origin(ν), μ_inner, x)
end

function _from_mvstd_with_rest_withorigin(::AbstractMeasure, ν_origin, μ_inner::StdMeasure, x)
    x_origin, x_rest = transport_from_mvstd_with_rest(ν_origin, x, μ_inner)
    from_origin(x_origin), x_rest
end

function _from_mvstd_with_rest_withorigin(ν::AbstractMeasure, NoTransportOrigin, μ_inner::StdMeasure, x)
    throw(ArgumentError("Don't know how to transport value of type $(nameof(typeof(x))) from power of $(nameof(typeof(μ_inner))) to $(nameof(typeof(ν)))"))
end


# Implement transport_to(NU::Type{<:StdMeasure}, μ) and transport_to(ν, MU::Type{<:StdMeasure})
# for user convenience:

_std_measure_for(::Type{M}, μ::Any) where {M<:StdMeasure} = _std_measure_for_impl(M, some_dof(μ))
_std_measure_for_impl(::Type{M}, ::StaticInteger{1}) where {M<:StdMeasure} = M()
_std_measure_for_impl(::Type{M}, dof::Integer) where {M<:StdMeasure} = M()^dof


function transport_to(ν, ::Type{MU}) where {MU<:StdMeasure}
    transport_to(ν, _std_measure_for(MU, ν))
end

function transport_to(::Type{NU}, μ) where {NU<:StdMeasure}
    transport_to(_std_measure_for(NU, μ), μ)
end


@inline transport_origin(μ::ProductMeasure) = _marginals_tp_origin(marginals(μ))
@inline from_origin(μ::ProductMeasure, x_origin) = _marginals_from_origin(marginals(μ), x_origin)

_marginals_tp_origin(::Ms) where Ms = NoTransportOrigin{PowerMeasure{M}}()


# Pull back from a product over a Fill to a power measure:

_marginals_tp_origin(marginals_μ::Fill) = marginals_μ.value^marginals_μ.axes
_marginals_from_origin(::Fill, x_origin) = x_origin


# Pull back from a NamedTuple product measure to a Tuple product measure:
#
# Maybe ToDo (breaking): For transport between NamedTuple-marginals we could
# match names where possible, even if given in different order, and transport
# between the remaining non-matching names in the order given. This may not
# be worth the additional complexity, though, since transport is typically
# used with a (power of a) standard measure on one side.

_marginals_tp_origin(marginals_μ::NamedTuple{names}) where names = productmeasure(values(marginals_μ))
_marginals_from_origin(::NamedTuple{names}, x_origin::NamedTuple) where names = _reorder_nt(x_origin, Val(names))


# Transport between two instances of ProductMeasure:

transport_def(ν::ProductMeasure, μ::ProductMeasure, x) = _marginal_transport_def(marginals(ν), marginals(μ), x)

function _marginal_transport_def(marginals_ν, marginals_μ, x)
    @assert size(marginals_ν) == size(marginals_μ) == size(x)  # Sanity check, should not fail
    transport_def.(marginals_ν, marginals_μ, x)
end

function _marginal_transport_def(marginals_ν::Tuple{Vararg{AbstractMeasure,N}}, marginals_μ::Tuple{Vararg{AbstractMeasure,N}}, x) where N
    @assert x isa Tuple{Vararg{AbstractMeasure,N}}  # Sanity check, should not fail
    map(transport_def, marginals_ν, marginals_μ, x)
end

function _marginal_transport_def(marginals_ν::AbstractVector{<:AbstractMeasure}, marginals_μ::Tuple{Vararg{AbstractMeasure,N}}, x) where N
    _marginal_transport_def(_as_tuple(marginals_ν, Val(N)), marginals_μ, x)
end

function _marginal_transport_def(marginals_ν::Tuple{Vararg{AbstractMeasure,N}}, marginals_μ::AbstractVector{<:AbstractMeasure}, x) where N
    _marginal_transport_def(marginals_ν, _as_tuple(marginals_μ, Val(N)), _as_tuple(x, Val(N)))
end



function transport_to_mvstd(ν_inner::StdMeasure, μ::ProductMeasure, ab)
    _marginals_to_mvstd(ν_inner, marginals(μ), ab)
end

function transport_from_mvstd_with_rest(ν::ProductMeasure, μ_inner::StdMeasure, x)
    a, x2 = transport_from_mvstd_with_rest(ν.α, μ_inner, x)
    b, x_rest = transport_from_mvstd_with_rest(ν.β, μ_inner, x2)
    return ν.f_c(a, b), x_rest
end


# Transport between a standard measure and Dirac:

@inline transport_from_mvstd_with_rest(ν::Dirac, ::StdMeasure, x::Any) = ν.x, x

@inline transport_to_mvstd(::StdMeasure, ::Dirac, ::Any) = Zeros{Bool}(0)



# Transport for products

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
