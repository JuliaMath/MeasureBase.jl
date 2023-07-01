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

function transport_def(ν::AbstractMeasure, μ::StdPowerMeasure{MU,1}, x) where MU
    μ_inner = _inner_stdmeasure(μ)
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
