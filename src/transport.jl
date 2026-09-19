"""
    f = transport_to(ν, μ)

Generates a [measurable function](https://en.wikipedia.org/wiki/Measurable_function)
`f` that transforms a value `x` distributed according to measure `μ` to
a value `y = f(x)` distributed according to a measure `ν`.

The [pushforward measure](https://en.wikipedia.org/wiki/Pushforward_measure)
from `μ` under `f` is equivalent to `ν`, so `f(rand(μ))` is equivalent
to `rand(ν)`. `f` supports `InverseFunctions.inverse` and
`ChangesOfVariables.with_logabsdet_jacobian`.

Measures are transported via standard measures: `x` is transported to the
standard measure type that the preferences of `ν` and `μ` promote to (see
[`MeasureBase.preferred_stdmeasure`](@ref)) and from there to `ν`.
Broadcasting `f` over an array of variates with flat storage, or over
the flat storage of a batch of variates (see
[`MeasureBase.mspace_ndims`](@ref)), transports the whole batch at once. A standard measure
type like `StdUniform` or `StdNormal` may also be used directly as the
source or target:

```julia
transport_to(StdNormal, μ)
transport_to(ν, StdNormal)
```

The transport partner is then an instance of the standard measure for
measures with scalar variates, and a power of it with as many elements as
the measure has degrees of freedom otherwise.

# Extended help

To support transport for a measure type, specialize
[`MeasureBase.transport_to_std`](@ref) and
[`MeasureBase.transport_from_std`](@ref) for its preferred standard measure
type, and declare [`MeasureBase.mspace_ndims`](@ref) for array variates.
Measure types with array variates should also implement the batched forms
[`MeasureBase.batched_transport_to_std`](@ref) and
[`MeasureBase.batched_transport_from_std`](@ref), which transport whole
batches of variates. Measures whose variates are composed of the variates
of other measures specialize the stream forms
[`MeasureBase.transport_to_std_with_rest`](@ref) and
[`MeasureBase.transport_from_std_with_rest`](@ref) instead (and their
batched forms). [`MeasureBase.transport_def`](@ref) may be specialized
for pairs of measure types with a direct transport.
"""
function transport_to end
export transport_to

"""
    transport_to(ν, μ, x)

Transport `x` from the measure `μ` to the measure `ν`, equivalent to
`transport_to(ν, μ)(x)`.

# Extended help

Variates of the right shape never throw: outside the support of `μ` the
result is `NaN` (elementwise for powers and products), values at the
boundary of the support may map to infinite values. Variates of the wrong
shape throw an `ArgumentError`. Transport implementations must not throw
outside the support, since the `NaN` masks evaluate both branches.
"""
transport_to(ν, μ, x) = transport_to(ν, μ)(x)

"""
    struct MeasureBase.TransportFunction <: Function

Transforms a variate from one measure to a variate of another.

In general `TransportFunction` should not be called directly, call
[`transport_to`](@ref) instead.
"""
struct TransportFunction{NU,MU} <: Function
    ν::NU
    μ::MU

    function TransportFunction{NU,MU}(ν::NU, μ::MU) where {NU,MU}
        return new{NU,MU}(ν, μ)
    end

    function TransportFunction(ν::NU, μ::MU) where {NU,MU}
        check_dof(ν, μ)
        return new{NU,MU}(ν, μ)
    end
end

@inline transport_to(ν, μ) = TransportFunction(asmeasure(ν), asmeasure(μ))

function Base.:(==)(a::TransportFunction, b::TransportFunction)
    return a.ν == b.ν && a.μ == b.μ
end

Base.@propagate_inbounds function (f::TransportFunction)(x)
    return transport_def(f.ν, f.μ, checked_arg(f.μ, x))
end

@inline function InverseFunctions.inverse(f::TransportFunction{NU,MU}) where {NU,MU}
    return TransportFunction{MU,NU}(f.μ, f.ν)
end

function ChangesOfVariables.with_logabsdet_jacobian(f::TransportFunction, x)
    y = f(x)
    logd_src = logdensityof(f.μ, x)
    logd_trg = logdensityof(f.ν, y)
    ladj = logd_src - logd_trg
    # Both densities being -Inf leaves the Jacobian undefined, zero is a safe choice then:
    fixed_ladj = ifelse(isneginf(logd_src) & isneginf(logd_trg), zero(ladj), ladj)
    return y, fixed_ladj
end

Base.:(∘)(::typeof(identity), f::TransportFunction) = f
Base.:(∘)(f::TransportFunction, ::typeof(identity)) = f

function Base.:∘(outer::TransportFunction, inner::TransportFunction)
    if !(outer.μ == inner.ν || isequal(outer.μ, inner.ν) || outer.μ ≈ inner.ν)
        throw(
            ArgumentError(
                "Cannot compose TransportFunction if source of outer doesn't equal target of inner.",
            ),
        )
    end
    return TransportFunction(outer.ν, inner.μ)
end

function Base.show(io::IO, f::TransportFunction)
    print(io, Base.typename(typeof(f)).name, "(")
    show(io, f.ν)
    print(io, ", ")
    show(io, f.μ)
    print(io, ")")
end

Base.show(io::IO, M::MIME"text/plain", f::TransportFunction) = show(io, f)


"""
    MeasureBase.transport_def(ν, μ, x)

Transport a variate `x` of `μ` to a variate of `ν`.

The default implementation transports `x` via the standard measure type the
preferences of `ν` and `μ` promote to. Specialize `transport_def` for pairs
of measure types with a direct transport.
"""
function transport_def end

@inline transport_def(ν, μ, x) = _transport_via_std(_transport_pivot(ν, μ), ν, μ, x)

function _transport_via_std(::Type{S}, ν, μ, x) where {S<:StdMeasure}
    z = transport_to_std(S, μ, x)
    y, z_rest = transport_from_std_with_rest(S, ν, _as_stdstream(z))
    if !isempty(z_rest)
        throw(ArgumentError("Degrees of freedom of source and target measure of a transport don't match"))
    end
    return y
end

@inline function _transport_pivot(ν, μ)
    _concrete_pivot(promote_stdmeasure(preferred_stdmeasure(ν), preferred_stdmeasure(μ)), ν, μ)
end
@inline function _concrete_pivot(::Type{S}, ν, μ) where {S<:StdMeasure}
    isconcretetype(S) || _throw_abstract_std(S)
    return S
end
@inline _concrete_pivot(::Type{AnyStdMeasure}, ν, μ) = StdUniform
function _concrete_pivot(::Type{<:NoStdTransport{MU}}, ν, μ) where {MU}
    throw(ArgumentError("No transport between measures of type $(nameof(typeof(ν))) and $(nameof(typeof(μ))), measures of type $(nameof(MU)) have no transport via standard measures"))
end

# Standard variates of scalar-variate measures are scalars, streams of
# standard variates are vectors:
@inline _as_stdstream(z::AbstractVector) = z
@inline _as_stdstream(z::Number) = SVector(z)


"""
    MeasureBase.transport_to_std(::Type{S}, μ, x)

Transport a variate `x` of `μ` to a variate of the standard measure type
`S`: a number if the variates of `μ` are scalars, a flat vector of length
[`getdof(μ)`](@ref) otherwise.

Measure types specialize `transport_to_std` for their preferred standard
measure type (see [`MeasureBase.preferred_stdmeasure`](@ref)), the generic
implementation converts between standard measure types.
"""
function transport_to_std end

@inline function transport_to_std(::Type{S}, μ, x) where {S<:StdMeasure}
    _to_std_via(S, preferred_stdmeasure(μ), μ, x)
end

@inline function _to_std_via(::Type{S}, ::Type{T}, μ, x) where {S<:StdMeasure,T<:StdMeasure}
    stdconvert(S, T, transport_to_std(T, μ, x))
end
function _to_std_via(::Type{S}, ::Type{S}, μ, x) where {S<:StdMeasure}
    throw(ArgumentError("Transport to $(nameof(S)) is not implemented for measures of type $(nameof(typeof(μ)))"))
end
function _to_std_via(::Type{S}, ::Type{AnyStdMeasure}, μ, x) where {S<:StdMeasure}
    throw(ArgumentError("Transport to standard measures is not implemented for measures of type $(nameof(typeof(μ)))"))
end
function _to_std_via(::Type{S}, ::Type, μ, x) where {S<:StdMeasure}
    throw(ArgumentError("Measures of type $(nameof(typeof(μ))) have no transport via standard measures"))
end


"""
    MeasureBase.transport_from_std(::Type{S}, μ, z)

Transport a variate `z` of the standard measure type `S` to a variate of
`μ`, the inverse of [`MeasureBase.transport_to_std`](@ref).
"""
function transport_from_std end

@inline function transport_from_std(::Type{S}, μ, z) where {S<:StdMeasure}
    _from_std_via(S, preferred_stdmeasure(μ), μ, z)
end

@inline function _from_std_via(::Type{S}, ::Type{T}, μ, z) where {S<:StdMeasure,T<:StdMeasure}
    transport_from_std(T, μ, stdconvert(T, S, z))
end
function _from_std_via(::Type{S}, ::Type{S}, μ, z) where {S<:StdMeasure}
    throw(ArgumentError("Transport from $(nameof(S)) is not implemented for measures of type $(nameof(typeof(μ)))"))
end
function _from_std_via(::Type{S}, ::Type{AnyStdMeasure}, μ, z) where {S<:StdMeasure}
    throw(ArgumentError("Transport from standard measures is not implemented for measures of type $(nameof(typeof(μ)))"))
end
function _from_std_via(::Type{S}, ::Type, μ, z) where {S<:StdMeasure}
    throw(ArgumentError("Measures of type $(nameof(typeof(μ))) have no transport via standard measures"))
end


"""
    MeasureBase.transport_to_std_with_rest(::Type{S}, μ, x)

Transport the variate of `μ` at the beginning of the stream `x` of
variate content to the standard measure type `S`.

Returns a tuple `(z, x_μ, x_rest)` of the flat vector `z` of standard
variates, the variate `x_μ` of `μ` consumed from the stream and the
unconsumed rest of the stream. See
[`MeasureBase.logdensityof_with_rest`](@ref) for the stream conventions.

The default implementation consumes a variate of the size given by
[`MeasureBase.mspace_flatsize`](@ref) or
[`MeasureBase.some_mspace_elsize`](@ref). Measure types whose variates are
composed of the variates of other measures implement
`transport_to_std_with_rest` instead of
[`MeasureBase.transport_to_std`](@ref).
"""
function transport_to_std_with_rest end

function transport_to_std_with_rest(::Type{S}, μ, x::AbstractVector) where {S<:StdMeasure}
    x_μ, x_rest = _consume_from_stream(x, _stream_consume_size(μ))
    return _as_stdstream(transport_to_std(S, μ, x_μ)), x_μ, x_rest
end

function transport_to_std_with_rest(::Type{S}, μ, x::NamedTuple) where {S<:StdMeasure}
    x_μ, x_rest = _split_after(x, Val(_mspace_names(μ)))
    return _as_stdstream(transport_to_std(S, μ, x_μ)), x_μ, x_rest
end


"""
    MeasureBase.transport_from_std_with_rest(::Type{S}, μ, z)

Transport the beginning of the flat stream `z` of standard variates of type
`S` to a variate of `μ`, consuming as many entries as `μ` requires.

Returns a tuple `(x, z_rest)` of the variate `x` and the unconsumed rest of
the stream. Measure types whose degrees of freedom depend on variate values
implement `transport_from_std_with_rest` instead of
[`MeasureBase.transport_from_std`](@ref).
"""
function transport_from_std_with_rest end

function transport_from_std_with_rest(::Type{S}, μ, z::AbstractVector) where {S<:StdMeasure}
    _from_std_with_rest_bydof(S, μ, z, fast_dof(μ))
end

function _from_std_with_rest_bydof(::Type{S}, μ, z::AbstractVector, n::IntegerLike) where {S}
    if maybestatic_length(z) < n
        throw(ArgumentError("Stream of standard variates too short during transport"))
    end
    z_μ, z_rest = _split_after(z, n)
    return transport_from_std(S, μ, _chunk_as_variate(μ, z_μ)), z_rest
end

function _from_std_with_rest_bydof(::Type{S}, μ, z::AbstractVector, ::AbstractNoDOF) where {S}
    throw(ArgumentError("Transport from standard measures requires measures of type $(nameof(typeof(μ))) to implement MeasureBase.transport_from_std_with_rest"))
end

# Scalar-variate measures take their standard variate as a number:
@inline _chunk_as_variate(μ, z) = _chunk_as_variate(z, _static_ndims(μ))
@inline _chunk_as_variate(z::AbstractVector, ::StaticInteger{0}) = z[begin]
@inline _chunk_as_variate(z::AbstractVector, ::Any) = z


"""
    transport_to(ν, ::Type{MU}) where {MU<:StdMeasure}
    transport_to(::Type{NU}, μ) where {NU<:StdMeasure}

As a user convenience, a standard measure type like [`StdUniform`](@ref),
[`StdExponential`](@ref), [`StdNormal`](@ref) or [`StdLogistic`](@ref)
may be used directly as the source or target of a measure transport.

The transport partner is an instance of the standard measure for measures
with scalar variates, and a power of it with
[`MeasureBase.some_dof(μ)`](@ref) (resp. `ν`) elements otherwise.
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

function _std_tp_partner(::Type{M}, μ) where {M<:StdMeasure}
    m = asmeasure(μ)
    _std_tp_partner_byrank(M, _static_ndims(m), m)
end
_std_tp_partner_byrank(::Type{M}, ::StaticInteger{0}, μ) where {M<:StdMeasure} = M()
_std_tp_partner_byrank(::Type{M}, ::Any, μ) where {M<:StdMeasure} = M()^some_dof(μ)


# Element-wise transport kernels for broadcasts and maps:
struct _ToStd{S} <: Function end
@inline (::_ToStd{S})(μ, x) where {S} = transport_to_std(S, μ, x)
struct _FromStd{S} <: Function end
@inline (::_FromStd{S})(μ, z) where {S} = transport_from_std(S, μ, z)

# Flat vector of standard variates from an array of standard variates:
@inline _flat_std_of(A::AbstractArray{<:Number}) = vec(A)
@inline _flat_std_of(A::AbstractArray{<:AbstractVector}) = _flatten_to_rv(vec(A))
@inline _flat_std_of(A::AbstractVector{<:AbstractVector}) = _flatten_to_rv(A)
