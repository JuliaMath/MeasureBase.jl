"""
    struct MeasureBase.NoTransportOrigin{NU}

Indicates that no (default) pullback measure is available for measures of
type `NU`.

See [`MeasureBase.transport_origin`](@ref).
"""
struct NoTransportOrigin{NU} end

"""
    MeasureBase.transport_origin(ν)

Default measure to pullback to resp. pushforward from when transforming
between `ν` and another measure.
"""
function transport_origin end

transport_origin(ν::NU) where {NU} = NoTransportOrigin{NU}()

"""
    MeasureBase.from_origin(ν, x)

Push `x` from `MeasureBase.transport_origin(μ)` forward to `ν`.
"""
function from_origin end

from_origin(ν::NU, ::Any) where {NU} = NoTransportOrigin{NU}()

"""
    MeasureBase.to_origin(ν, y)

Pull `y` from `ν` back to `MeasureBase.transport_origin(ν)`.
"""
function to_origin end

to_origin(ν::NU, ::Any) where {NU} = NoTransportOrigin{NU}()

"""
    struct MeasureBase.NoTransport{NU,MU} end

Indicates that no transformation from a measure of type `MU` to a measure of
type `NU` could be found.
"""
struct NoTransport{NU,MU} end

"""
    f = transport_to(ν, μ)

Generates a [measurable function](https://en.wikipedia.org/wiki/Measurable_function)
`f` that transforms a value `x` distributed according to measure `μ` to
a value `y = f(x)` distributed according to a measure `ν`.

The [pushforward measure](https://en.wikipedia.org/wiki/Pushforward_measure)
from `μ` under `f` is is equivalent to `ν`.

If terms of random values this implies that `f(rand(μ))` is equivalent to
`rand(ν)` (if `rand(μ)` and `rand(ν)` are supported).

The resulting function `f` should support
`ChangesOfVariables.with_logabsdet_jacobian(f, x)` if mathematically well-defined,
so that densities of `ν` can be derived from densities of `μ` via `f` (using
appropriate base measures).

Returns NoTransportOrigin{typeof(ν),typeof(μ)} if no transformation from
`μ` to `ν` can be found.

To add transformation rules for a measure type `MyMeasure`, specialize

* `MeasureBase.transport_def(ν::SomeStdMeasure, μ::CustomMeasure, x) = ...`
* `MeasureBase.transport_def(ν::MyMeasure, μ::SomeStdMeasure, x) = ...`

and/or

* `MeasureBase.transport_origin(ν::MyMeasure) = SomeMeasure(...)`
* `MeasureBase.from_origin(μ::MyMeasure, x) = y`
* `MeasureBase.to_origin(μ::MyMeasure, y) = x`

and ensure `MeasureBase.getdof(μ::MyMeasure)` is defined correctly.

A standard measure type like `StdUniform`, `StdExponential` or
`StdLogistic` may also be used as the source or target of the transform:

```julia
f_to_uniform(StdUniform, μ)
f_to_uniform(ν, StdUniform)
```

Depending on [`getdof(μ)`](@ref) (resp. `ν`), an instance of the standard
distribution itself or a power of it (e.g. `StdUniform()` or
`StdUniform()^dof`) will be chosen as the transformation partner.
"""
function transport_to end

@inline transport_to(ν, μ) = TransportFunction(asmeasure(ν), asmeasure(μ))

"""
    transport_to(ν, μ, x)

Transport `x` from the measure `μ` to the measure `ν`
"""
transport_to(ν, μ, x) = transport_to(ν, μ)(x)


"""
    transport_def(ν, μ, x)

Transforms a value `x` distributed according to `μ` to a value `y` distributed
according to `ν`.

If no specialized `transport_def(::MU, ::NU, ...)` is available then
the default implementation of`transport_def(ν, μ, x)` uses the following
strategy:

* Evaluate [`transport_origin`](@ref) for μ and ν. Transform between
  each and it's origin, if available, and use the origin(s) as intermediate
  measures for another transformation.

* If all else fails, try to transform from μ to a standard multivariate
  uniform measure and then to ν.

See [`transport_to`](@ref).
"""
function transport_def end

function transport_def(ν, μ, x)
    _transport_between_origins(ν, _origin_depth(ν), _origin_depth(μ), μ, x)
end

@inline function _origin_depth(ν::NU) where {NU}
    ν_0 = ν
    Base.Cartesian.@nexprs 10 i -> begin  # 10 is just some "big enough" number
        ν_{i} = transport_origin(ν_{i - 1})
        if ν_{i} isa NoTransportOrigin
            return static(i - 1)
        end
    end
    return static(10)
end

_origin_depth_pullback(ΔΩ) = NoTangent(), NoTangent()
ChainRulesCore.rrule(::typeof(_origin_depth), ν) = _origin_depth(ν), _origin_depth_pullback

# If both both measures have no origin:
function _transport_between_origins(ν, ::StaticInteger{0}, ::StaticInteger{0}, μ, x)
    _transport_with_intermediate(ν, _transport_intermediate(ν, μ), μ, x)
end

@generated function _transport_between_origins(
    ν,
    ::StaticInteger{n_ν},
    ::StaticInteger{n_μ},
    μ,
    x,
) where {n_ν,n_μ}
    prog = quote
        μ0 = μ
        x0 = x
        ν0 = ν
    end
    for i in 1:n_μ
        μ_i = Symbol(:μ, i)
        μ_last = Symbol(:μ, i - 1)
        push!(prog.args, :($μ_i = transport_origin($μ_last)))
    end
    for i in 1:n_μ
        x_i = Symbol(:x, i)
        x_last = Symbol(:x, i - 1)
        μ_last = Symbol(:μ, i - 1)
        push!(prog.args, :($x_i = to_origin($μ_last, $x_last)))
    end
    for i in 1:(n_ν)
        ν_i = Symbol(:ν, i)
        ν_last = Symbol(:ν, i - 1)
        push!(prog.args, :($ν_i = transport_origin($ν_last)))
    end
    μ_im = Symbol(:μ, n_μ)
    x_im = Symbol(:x, n_μ)
    ν_im = Symbol(:ν, n_ν)
    y_im = Symbol(:y, n_ν)
    push!(prog.args, :($y_im = transport_def($ν_im, $μ_im, $x_im)))
    for i in (n_ν-1):-1:0
        y_i = Symbol(:y, i)
        y_last = Symbol(:y, i + 1)
        ν_last = Symbol(:ν, i)
        push!(prog.args, :($y_i = from_origin($ν_last, $y_last)))
    end
    push!(prog.args, :(return y0))
    return prog
end

@inline _transport_intermediate(ν, μ) = _transport_intermediate(getdof(ν), getdof(μ))
@inline _transport_intermediate(::Integer, n_μ::Integer) = StdUniform()^n_μ
@inline _transport_intermediate(::StaticInteger{1}, ::StaticInteger{1}) = StdUniform()

_call_transport_def(ν, μ, x) = transport_def(ν, μ, x)
_call_transport_def(::Any, ::Any, x::NoTransportOrigin) = x
_call_transport_def(::Any, ::Any, x::NoTransport) = x

function _transport_with_intermediate(ν, m, μ, x)
    z = _call_transport_def(m, μ, x)
    y = _call_transport_def(ν, m, z)
    return y
end

# Prevent infinite recursion in case vartransform_intermediate doesn't change type:
@inline function _transport_with_intermediate(::NU, ::NU, ::MU, ::Any) where {NU,MU}
    NoTransport{NU,MU}()
end
@inline function _transport_with_intermediate(::NU, ::MU, ::MU, ::Any) where {NU,MU}
    NoTransport{NU,MU}()
end

"""
    struct TransportFunction <: Function

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

function Base.:(==)(a::TransportFunction, b::TransportFunction)
    return a.ν == b.ν && a.μ == b.μ
end

Base.@propagate_inbounds function (f::TransportFunction)(x)
    return _call_transport_def(f.ν, f.μ, checked_arg(f.μ, x))
end

@inline function InverseFunctions.inverse(f::TransportFunction{NU,MU}) where {NU,MU}
    return TransportFunction{MU,NU}(f.μ, f.ν)
end

function ChangesOfVariables.with_logabsdet_jacobian(f::TransportFunction, x)
    y = f(x)
    logpdf_src = logdensityof(f.μ, x)
    logpdf_trg = logdensityof(f.ν, y)
    ladj = logpdf_src - logpdf_trg
    # If logpdf_src and logpdf_trg are -Inf setting lafj to zero is safe:
    fixed_ladj = logpdf_src == logpdf_trg == -Inf ? zero(ladj) : ladj
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
    abstract type TransformVolCorr

Provides control over density correction by transform volume element.
Either [`NoVolCorr()`](@ref) or [`WithVolCorr()`](@ref)
"""
abstract type TransformVolCorr end

"""
    NoVolCorr()

Indicate that density calculations should ignore the volume element of
variate transformations. Should only be used in special cases in which
the volume element has already been taken into account in a different
way.
"""
struct NoVolCorr <: TransformVolCorr end

"""
    WithVolCorr()

Indicate that density calculations should take the volume element of
variate transformations into account (typically via the
log-abs-det-Jacobian of the transform).
"""
struct WithVolCorr <: TransformVolCorr end
