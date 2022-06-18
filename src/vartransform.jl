"""
    struct MeasureBase.NoTransformOrigin{NU}

Indicates that no (default) pullback measure is available for measures of
type `NU`.

See [`MeasureBase.vartransform_origin`](@ref).
"""
struct NoTransformOrigin{NU} end


"""
    MeasureBase.vartransform_origin(ν)

Default measure to pullback to resp. pushforward from when transforming
between `ν` and another measure.
"""
function vartransform_origin end

vartransform_origin(ν::NU) where NU = NoTransformOrigin{NU}()


"""
    MeasureBase.from_origin(ν, x)

Push `x` from `MeasureBase.vartransform_origin(μ)` forward to `ν`.
"""
function from_origin end

from_origin(ν::NU, ::Any) where NU = NoTransformOrigin{NU}()


"""
    MeasureBase.to_origin(ν, y)

Pull `y` from `ν` back to `MeasureBase.vartransform_origin(ν)`.
"""
function to_origin end

to_origin(ν::NU, ::Any) where NU = NoTransformOrigin{NU}(ν)


"""
    struct MeasureBase.NoVarTransform{NU,MU} end

Indicates that no transformation from a measure of type `MU` to a measure of
type `NU` could be found.
"""
struct NoVarTransform{NU,MU} end


"""
    f = vartransform(ν, μ)

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

Returns NoTransformOrigin{typeof(ν),typeof(μ)} if no transformation from
`μ` to `ν` can be found.

To add transformation rules for a measure type `MyMeasure`, specialize

* `MeasureBase.vartransform_def(ν::SomeStdMeasure, μ::CustomMeasure, x) = ...`
* `MeasureBase.vartransform_def(ν::MyMeasure, μ::SomeStdMeasure, x) = ...`

and/or

* `MeasureBase.vartransform_origin(ν::MyMeasure) = SomeMeasure(...)`
* `MeasureBase.from_origin(μ::MyMeasure, y) = x`
* `MeasureBase.to_origin(μ::MyMeasure, x) = y`

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
function vartransform end


"""
    vartransform_def(ν, μ, x)

Transforms a value `x` distributed according to `μ` to a value `y` distributed
according to `ν`.

If no specialized `vartransform_def(::MU, ::NU, ...)` is available then
the default implementation of`vartransform_def(ν, μ, x)` uses the following
strategy:

* Evaluate [`vartransform_origin`](@ref) for μ and ν. Transform between
  each and it's origin, if available, and use the origin(s) as intermediate
  measures for another transformation.

* If all else fails, try to transform from μ to a standard multivariate
  uniform measure and then to ν.

See [`vartransform`](@ref).
"""
function vartransform_def end

vartransform_def(::Any, ::Any, x::NoTransformOrigin) = x
vartransform_def(::Any, ::Any, x::NoVarTransform) = x

function vartransform_def(ν, μ, x)
    check_dof(ν, μ)
    _vartransform_with_intermediate(ν, _checked_vartransform_origin(ν), _checked_vartransform_origin(μ), μ, x)
end


@inline _origin_must_have_separate_type(::Type{MU}, μ_o) where MU = μ_o
function _origin_must_have_separate_type(::Type{MU}, μ_o::MU) where MU
    throw(ArgumentError("Measure of type $MU and its origin must have separate types"))
end

@inline function _checked_vartransform_origin(μ::MU) where MU
    μ_o = vartransform_origin(μ)
    _origin_must_have_separate_type(MU, μ_o)
end


function _vartransform_with_intermediate(ν, ν_o, μ_o, μ, x)
    x_o = to_origin(μ, x)
    # If μ is a pushforward then checked_var may have been bypassed, so check now:
    y_o = vartransform_def(ν_o, μ_o, checked_var(μ_o, x_o))
    y = from_origin(ν, y_o)
    return y
end

function _vartransform_with_intermediate(ν, ν_o, ::NoTransformOrigin, μ, x)
    y_o = vartransform_def(ν_o, μ, x)
    y = from_origin(ν, y_o)
    return y
end

function _vartransform_with_intermediate(ν, ::NoTransformOrigin, μ_o, μ, x)
    x_o = to_origin(μ, x)
    # If μ is a pushforward then checked_var may have been bypassed, so check now:
    y = vartransform_def(ν, μ_o, checked_var(μ_o, x_o))
    return y
end

function _vartransform_with_intermediate(ν, ::NoTransformOrigin, ::NoTransformOrigin, μ, x)
    _vartransform_with_intermediate(ν, _vartransform_intermediate(ν, μ), μ, x)
end


@inline _vartransform_intermediate(ν, μ) = _vartransform_intermediate(getdof(ν), getdof(μ))
@inline _vartransform_intermediate(::Integer, n_μ::Integer) = StdUniform()^n_μ
@inline _vartransform_intermediate(::StaticInt{1}, ::StaticInt{1}) = StdUniform()

function _vartransform_with_intermediate(ν, m, μ, x)
    z = vartransform_def(m, μ, x)
    y = vartransform_def(ν, m, z)
    return y
end

# Prevent infinite recursion in case vartransform_intermediate doesn't change type:
@inline _vartransform_with_intermediate(::NU, ::NU, ::MU, ::Any) where {NU,MU} = NoVarTransform{NU,MU}()
@inline _vartransform_with_intermediate(::NU, ::MU, ::MU, ::Any) where {NU,MU} = NoVarTransform{NU,MU}()


"""
    struct VarTransformation <: Function

Transforms a variate from one measure to a variate of another.

In general `VarTransformation` should not be called directly, call
[`vartransform`](@ref) instead.
"""
struct VarTransformation{NU,MU} <: Function
    ν::NU
    μ::MU

    function VarTransformation{NU,MU}(ν::NU, μ::MU) where {NU,MU}
        return new{NU,MU}(ν, μ)
    end

    function VarTransformation(ν::NU, μ::MU) where {NU,MU}
        check_dof(ν, μ)
        return new{NU,MU}(ν, μ)
    end
end

@inline vartransform(ν, μ) = VarTransformation(ν, μ)

function Base.:(==)(a::VarTransformation, b::VarTransformation)
    return a.ν == b.ν && a.μ == b.μ
end


Base.@propagate_inbounds function (f::VarTransformation)(x)
    return vartransform_def(f.ν, f.μ, checked_var(f.μ, x))
end

@inline function InverseFunctions.inverse(f::VarTransformation{NU,MU}) where {NU,MU}
    return VarTransformation{MU,NU}(f.μ, f.ν)
end


function ChangesOfVariables.with_logabsdet_jacobian(f::VarTransformation, x)
    y = f(x)
    logpdf_src = logdensityof(f.μ, x)
    logpdf_trg = logdensityof(f.ν, y)
    ladj = logpdf_src - logpdf_trg
    # If logpdf_src and logpdf_trg are -Inf setting lafj to zero is safe:
    fixed_ladj = logpdf_src == logpdf_trg == -Inf ? zero(ladj) : ladj
    return y, fixed_ladj
end


Base.:(∘)(::typeof(identity), f::VarTransformation) = f
Base.:(∘)(f::VarTransformation, ::typeof(identity)) = f

function Base.:∘(outer::VarTransformation, inner::VarTransformation)
    if !(outer.μ == inner.ν || isequal(outer.μ, inner.ν) || outer.μ ≈ inner.ν)
        throw(ArgumentError("Cannot compose VarTransformation if source of outer doesn't equal target of inner."))
    end 
    return VarTransformation(outer.ν, inner.μ)
end


function Base.show(io::IO, f::VarTransformation)
    print(io, Base.typename(typeof(f)).name, "(")
    show(io, f.ν)
    print(io, ", ")
    show(io, f.μ)
    print(io, ")")
end

Base.show(io::IO, M::MIME"text/plain", f::VarTransformation) = show(io, f)
