"""
    transport_to(ν, ::Type{MU}) where {MU<:StdMeasure}
    transport_to(::Type{NU}, μ) where {NU<:StdMeasure}

As a user convenience, a standard measure type like [`StdUniform`](@ref),
[`StdExponential`](@ref), [`StdNormal`](@ref) or [`StdLogistic`](@ref)
may be used directly as the source or target of a measure transport.

Depending on [`MeasureBase.some_dof(μ)`](@ref) (resp. `ν`), an instance of
the standard measure itself or a power of it will be automatically chosen as
the transport partner.

Example:

```julia
transport_to(StdNormal, μ)
transport_to(ν, StdNormal)
```
"""
function transport_to(ν, ::Type{MU}) where {MU<:StdMeasure}
    transport_to(ν, _std_tp_partner(MU, ν))
end

function transport_to(::Type{NU}, μ) where {NU<:StdMeasure}
    transport_to(_std_tp_partner(NU, μ), μ)
end

function transport_to(::Type{NU}, ::Type{MU}) where {NU<:StdMeasure,MU<:StdMeasure}
    throw(
        ArgumentError(
            "Can't construct a transport function between the types of two standard measures, need a measure instance on one side",
        ),
    )
end

_std_tp_partner(::Type{M}, μ) where {M<:StdMeasure} = _std_tp_partner_bydof(M, some_dof(μ))
_std_tp_partner_bydof(::Type{M}, ::StaticInteger{1}) where {M<:StdMeasure} = M()
_std_tp_partner_bydof(::Type{M}, dof::IntegerLike) where {M<:StdMeasure} = M()^dof
function _std_tp_partner_bydof(::Type{M}, ::AbstractNoDOF{MU}) where {M<:StdMeasure,MU}
    throw(
        ArgumentError(
            "Can't determine a standard transport partner for measures of type $(nameof(MU))",
        ),
    )
end


# For transport, always pull a multi-dimensional PowerMeasure back to a
# one-dimensional PowerMeasure first:

const _PowerMeasureRank1{M} = PowerMeasure{M,<:NTuple{1,OneToLike}}

function transport_origin(μ::PowerMeasure)
    pwr_base(μ)^prod(pwr_size(μ))
end

function to_origin(μ::PowerMeasure, x)
    maybestatic_reshape(x, (prod(pwr_size(μ)),))
end

function from_origin(μ::PowerMeasure, x_origin)
    # Sanity check, should never fail:
    @assert x_origin isa AbstractVector
    maybestatic_reshape(x_origin, pwr_size(μ))
end


# A one-dimensional PowerMeasure has an origin if its parent has an origin:

function transport_origin(μ::_PowerMeasureRank1)
    _pwr_origin(typeof(μ), transport_origin(pwr_base(μ)), pwr_axes(μ))
end
_pwr_origin(::Type{MU}, parent_origin, axes) where {MU} = parent_origin^axes
_pwr_origin(::Type{MU}, ::NoTransportOrigin, axes) where {MU} = NoTransportOrigin{MU}()

function to_origin(μ::_PowerMeasureRank1, x)
    to_origin.(Ref(pwr_base(μ)), x)
end

function from_origin(μ::_PowerMeasureRank1, x_origin)
    # Sanity check, should never fail:
    @assert x_origin isa AbstractVector
    from_origin.(Ref(pwr_base(μ)), x_origin)
end


# Transport between powers of standard measures, of any rank:

function _stdpow_transport_def(ν::StdPowerMeasure{NU}, μ::StdPowerMeasure{MU}, x) where {NU,MU}
    y = transport_to(pwr_base(ν), pwr_base(μ)).(x)
    maybestatic_reshape(y, pwr_size(ν))
end

function _stdpow_transport_def(ν::StdPowerMeasure{MU}, μ::StdPowerMeasure{MU}, x) where {MU}
    maybestatic_reshape(x, pwr_size(ν))
end

transport_def(ν::StdPowerMeasure{NU}, μ::StdPowerMeasure{MU}, x) where {NU,MU} =
    _stdpow_transport_def(ν, μ, x)

# Disambiguation with the mvstd gateway methods below:
transport_def(ν::StdPowerMeasure{NU,1}, μ::StdPowerMeasure{MU}, x) where {NU,MU} =
    _stdpow_transport_def(ν, μ, x)
transport_def(ν::StdPowerMeasure{NU}, μ::StdPowerMeasure{MU,1}, x) where {NU,MU} =
    _stdpow_transport_def(ν, μ, x)
transport_def(ν::StdPowerMeasure{NU,1}, μ::StdPowerMeasure{MU,1}, x) where {NU,MU} =
    _stdpow_transport_def(ν, μ, x)


# Transport between univariate standard measures and one-dimensional power
# measures of size one:

function transport_def(ν::StdMeasure, μ::StdPowerMeasure{MU,1}, x) where {MU}
    return transport_def(ν, pwr_base(μ), only(x))
end

function transport_def(ν::StdPowerMeasure{NU,1}, μ::StdMeasure, x) where {NU}
    sz_ν = pwr_size(ν)
    @assert prod(sz_ν) == 1
    return maybestatic_fill(transport_def(pwr_base(ν), μ, x), sz_ν)
end


# Transport to a multivariate standard measure from any measure:

function transport_def(ν::StdPowerMeasure{NU,1}, μ::AbstractMeasure, x) where {NU}
    y = transport_to_mvstd(pwr_base(ν), μ, x)
    if maybestatic_length(y) != maybestatic_length(ν)
        throw(ArgumentError("Length of transport target doesn't match variate DOF during transport"))
    end
    return y
end

function transport_to_mvstd(ν_inner::StdMeasure, μ::AbstractMeasure, x)
    return _to_mvstd_withdof(ν_inner, μ, fast_dof(μ), x)
end

# For standard measures and their powers specialized `transport_def` methods
# exist, for other measures the origin-based machinery must be used directly
# instead of `transport_def`, to prevent infinite dispatch recursion via the
# gateway methods above:
const _StdOrStdPowerMeasure = Union{StdMeasure,StdPowerMeasure}

_transport_def_nongateway(ν::_StdOrStdPowerMeasure, μ::_StdOrStdPowerMeasure, x) =
    transport_def(ν, μ, x)
function _transport_def_nongateway(ν, μ, x)
    _transport_between_origins(ν, _origin_depth(ν), _origin_depth(μ), μ, x)
end

function _to_mvstd_withdof(ν_inner::StdMeasure, μ::AbstractMeasure, dof_μ::IntegerLike, x)
    _transport_def_nongateway(ν_inner^dof_μ, μ, x)
end

function _to_mvstd_withdof(ν_inner::StdMeasure, μ::AbstractMeasure, ::AbstractNoDOF, x)
    _to_mvstd_withorigin(ν_inner, μ, transport_origin(μ), x)
end

function _to_mvstd_withorigin(ν_inner::StdMeasure, μ::AbstractMeasure, μ_origin, x)
    x_origin = to_origin(μ, x)
    transport_to_mvstd(ν_inner, μ_origin, x_origin)
end

function _to_mvstd_withorigin(ν_inner::StdMeasure, μ::AbstractMeasure, ::NoTransportOrigin, x)
    throw(
        ArgumentError(
            "Don't know how to transport values of type $(nameof(typeof(x))) from $(nameof(typeof(μ))) to a power of $(nameof(typeof(ν_inner)))",
        ),
    )
end


# Transport from a multivariate standard measure to any measure:

function transport_def(ν::AbstractMeasure, μ::StdPowerMeasure{MU,1}, x) where {MU}
    _transport_from_mvstd(ν, pwr_base(μ), x)
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
    return _from_mvstd_with_rest_withdof(ν, dof_ν, μ_inner, x)
end

function _from_mvstd_with_rest_withdof(
    ν::AbstractMeasure,
    dof_ν::IntegerLike,
    μ_inner::StdMeasure,
    x,
)
    len_x = maybestatic_length(x)

    # Since we can't check the DOF of measures like Bind, we could "run out
    # of x" if the original x was too short. `transport_to` below will detect
    # this, but better throw a more informative exception here:
    if len_x < dof_ν
        throw(ArgumentError("Variate too short during transport"))
    end

    x_inner_dof, x_rest = _split_after(x, dof_ν)
    y = _transport_def_nongateway(ν, μ_inner^dof_ν, x_inner_dof)
    return y, x_rest
end

function _from_mvstd_with_rest_withdof(
    ν::AbstractMeasure,
    ::AbstractNoDOF,
    μ_inner::StdMeasure,
    x,
)
    _from_mvstd_with_rest_withorigin(ν, transport_origin(ν), μ_inner, x)
end

function _from_mvstd_with_rest_withorigin(
    ν::AbstractMeasure,
    ν_origin,
    μ_inner::StdMeasure,
    x,
)
    x_origin, x_rest = transport_from_mvstd_with_rest(ν_origin, μ_inner, x)
    from_origin(ν, x_origin), x_rest
end

function _from_mvstd_with_rest_withorigin(
    ν::AbstractMeasure,
    ::NoTransportOrigin,
    μ_inner::StdMeasure,
    x,
)
    throw(
        ArgumentError(
            "Don't know how to transport a value of type $(nameof(typeof(x))) from a power of $(nameof(typeof(μ_inner))) to $(nameof(typeof(ν)))",
        ),
    )
end


# Transport between a standard measure and Dirac:

@inline transport_from_mvstd_with_rest(ν::Dirac, ::StdMeasure, x::Any) = ν.x, x

@inline transport_to_mvstd(::StdMeasure, ::Dirac, ::Any) = FillArrays.Zeros{Bool}(0)


# Pull back from a product over a Fill to a power measure:

@inline transport_origin(μ::ProductMeasure) = _marginals_tp_origin(marginals(μ))
@inline to_origin(μ::ProductMeasure, x) = _marginals_to_origin(marginals(μ), x)
@inline from_origin(μ::ProductMeasure, x_origin) =
    _marginals_from_origin(marginals(μ), x_origin)

_marginals_tp_origin(::Ms) where {Ms} = NoTransportOrigin{ProductMeasure{Ms}}()

_marginals_tp_origin(marginals_μ::FillArrays.Fill) =
    _fill_value(marginals_μ)^_fill_axes(marginals_μ)
_marginals_to_origin(::FillArrays.Fill, x) = x
_marginals_from_origin(::FillArrays.Fill, x_origin) = x_origin


# Pull back from a NamedTuple product measure to a Tuple product measure:
#
# Maybe ToDo (breaking): For transport between NamedTuple-marginals we could
# match names where possible, even if given in different order, and transport
# between the remaining non-matching names in the order given. This may not
# be worth the additional complexity, though, since transport is typically
# used with a (power of a) standard measure on one side.

_marginals_tp_origin(marginals_μ::NamedTuple{names}) where {names} =
    productmeasure(values(marginals_μ))
_marginals_to_origin(::NamedTuple{names}, x::NamedTuple{names}) where {names} = values(x)
_marginals_from_origin(::NamedTuple{names}, x_origin::Tuple) where {names} =
    NamedTuple{names}(x_origin)


# Transport between two instances of ProductMeasure:

transport_def(ν::ProductMeasure, μ::ProductMeasure, x) =
    _marginal_transport_def(marginals(ν), marginals(μ), x)

function _marginal_transport_def(marginals_ν, marginals_μ, x)
    @assert size(marginals_ν) == size(marginals_μ) == size(x)  # Sanity check, should not fail
    transport_def.(marginals_ν, marginals_μ, x)
end

function _marginal_transport_def(
    marginals_ν::Tuple{Vararg{AbstractMeasure,N}},
    marginals_μ::Tuple{Vararg{AbstractMeasure,N}},
    x,
) where {N}
    map(transport_def, marginals_ν, marginals_μ, x)
end

function _marginal_transport_def(
    marginals_ν::AbstractVector{<:AbstractMeasure},
    marginals_μ::Tuple{Vararg{AbstractMeasure,N}},
    x,
) where {N}
    _marginal_transport_def(_as_tuple(marginals_ν, Val(N)), marginals_μ, x)
end

function _marginal_transport_def(
    marginals_ν::Tuple{Vararg{AbstractMeasure,N}},
    marginals_μ::AbstractVector{<:AbstractMeasure},
    x,
) where {N}
    _marginal_transport_def(marginals_ν, _as_tuple(marginals_μ, Val(N)), _as_tuple(x, Val(N)))
end


# Transport from a ProductMeasure to a standard measure:

function transport_to_mvstd(ν_inner::StdMeasure, μ::ProductMeasure, x)
    _marginals_to_mvstd(ν_inner, marginals(μ), x)
end

struct _TransportToMvStd{NU<:StdMeasure} <: Function end
(::_TransportToMvStd{NU})(μ, x) where {NU} = transport_to_mvstd(NU(), μ, x)

function _marginals_to_mvstd(::NU, marginals_μ::Tuple, x::Tuple) where {NU<:StdMeasure}
    _flatten_to_rv(map(_TransportToMvStd{NU}(), marginals_μ, x))
end

function _marginals_to_mvstd(
    ν::NU,
    marginals_μ::NamedTuple{names},
    x::NamedTuple{names},
) where {NU<:StdMeasure,names}
    _marginals_to_mvstd(ν, values(marginals_μ), values(x))
end

function _marginals_to_mvstd(::NU, marginals_μ, x) where {NU<:StdMeasure}
    _flatten_to_rv(broadcast(_TransportToMvStd{NU}(), marginals_μ, x))
end


# Transport from a standard measure to a ProductMeasure, with rest:

const _MaybeUnknownDOF = Union{IntegerLike,AbstractNoDOF}

const _KnownDOFs = Union{Tuple{Vararg{IntegerLike,N}} where N,StaticVector{<:IntegerLike}}

function transport_from_mvstd_with_rest(ν::ProductMeasure, μ_inner::StdMeasure, x)
    νs = marginals(ν)
    dofs = map(fast_dof, νs)
    return _marginals_from_mvstd_with_rest(νs, dofs, μ_inner, x)
end

function transport_from_mvstd_with_rest(
    ν::ProductMeasure{<:NamedTuple{names}},
    μ_inner::StdMeasure,
    x,
) where {names}
    ys, x_rest =
        transport_from_mvstd_with_rest(productmeasure(values(marginals(ν))), μ_inner, x)
    return NamedTuple{names}(ys), x_rest
end

function _dof_access_firstidxs(dofs::Tuple{Vararg{IntegerLike,N}}, first_idx) where {N}
    cumsum((first_idx, dofs[begin:(end-1)]...))
end

function _dof_access_firstidxs(dofs::AbstractVector{<:IntegerLike}, first_idx)
    # ToDo: Improve implementation (reduce memory allocations):
    cumsum(vcat([eltype(dofs)(first_idx)], dofs[begin:(end-1)]))
end

function _split_x_by_marginals_with_rest(
    dofs::Union{Tuple,AbstractVector},
    x::AbstractVector{<:Real},
)
    x_idxs = maybestatic_eachindex(x)
    first_idxs = _dof_access_firstidxs(dofs, maybestatic_first(x_idxs))
    xs = map((from, n) -> _get_or_view(x, from, from + n - one(n)), first_idxs, dofs)
    x_rest = _get_or_view(x, first_idxs[end] + dofs[end], maybestatic_last(x_idxs))
    return xs, x_rest
end

function _marginals_from_mvstd_with_rest(
    νs,
    dofs::_KnownDOFs,
    μ_inner::StdMeasure,
    x::AbstractVector{<:Real},
)
    xs, x_rest = _split_x_by_marginals_with_rest(dofs, x)
    μs = map(n -> μ_inner^n, dofs)
    ys = map(transport_def, νs, μs, xs)
    return ys, x_rest
end

function _marginals_from_mvstd_with_rest(
    νs,
    dofs,
    μ_inner::StdMeasure,
    x::AbstractVector{<:Real},
)
    _marginals_from_mvstd_with_rest_nodof(νs, μ_inner, x)
end

function _marginals_from_mvstd_with_rest_nodof(
    νs::Tuple{Vararg{AbstractMeasure}},
    μ_inner::StdMeasure,
    x::AbstractVector{<:Real},
)
    # ToDo: Check for type stability, may need a generated function:
    y1, x_rest = transport_from_mvstd_with_rest(νs[1], μ_inner, x)
    y2_end, x_final_rest = _marginals_from_mvstd_with_rest_nodof(Base.tail(νs), μ_inner, x_rest)
    return (y1, y2_end...), x_final_rest
end

_marginals_from_mvstd_with_rest_nodof(::Tuple{}, ::StdMeasure, x::AbstractVector{<:Real}) =
    (), x

function _marginals_from_mvstd_with_rest_nodof(
    νs::AbstractVector{M},
    μ_inner::StdMeasure,
    x::AbstractVector{<:Real},
) where {M<:AbstractMeasure}
    if isconcretetype(M)
        # Marginals of concrete type produce variates of uniform type, so
        # the loop below is type stable (the type of the remaining stream
        # stays invariant under repeated view-taking):
        idxs = eachindex(νs)
        y1, x_rest = transport_from_mvstd_with_rest(νs[first(idxs)], μ_inner, x)
        ys = Vector{typeof(y1)}(undef, length(idxs))
        ys[begin] = y1
        j = firstindex(ys) + 1
        for i in Iterators.drop(idxs, 1)
            ys[j], x_rest = transport_from_mvstd_with_rest(νs[i], μ_inner, x_rest)
            j += 1
        end
        return ys, x_rest
    else
        # Fallback for marginals of mixed type:
        ys_any = Vector{Any}(undef, length(eachindex(νs)))
        x_rest = x
        for (i, ν) in zip(eachindex(ys_any), νs)
            ys_any[i], x_rest = transport_from_mvstd_with_rest(ν, μ_inner, x_rest)
        end
        return [y for y in ys_any], x_rest
    end
end
