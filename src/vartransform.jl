"""
    struct MeasureBase.NoTransformOrigin{MU}

Indicates that no (default) pullback measure is available for measures of
type `MU`.

See [`MeasureBase.vartransform_origin`](@ref).
"""
struct NoTransformOrigin{MU} end


"""
    MeasureBase.vartransform_origin(μ)

Default measure to pullback to resp. pushforward from when transforming
between `μ` and another measure.
"""
function vartransform_origin end

vartransform_origin(m::M) where M = NoTransformOrigin{M}()


"""
    MeasureBase.from_origin(μ, y)

Push `y` from `MeasureBase.vartransform_origin(μ)` forward to `μ`.
"""
function from_origin end

from_origin(m::M) where M = NoTransformOrigin{M}()


"""
    MeasureBase.to_origin(μ, x)

Pull `x` from `μ` back to `MeasureBase.vartransform_origin(μ)`.
"""
function to_origin end

to_origin(m::M) where M = NoTransformOrigin{M}()


"""
    struct MeasureBase.NoVarTransform{NU,MU} end

Indicates that no transformation from a measure of type `MU` to a measure of
type `NU` could be found.
"""
struct NoVarTransform{NU,MU} end


"""
    f = vartransform(ν, μ)

Generates a [measurable function](https://en.wikipedia.org/wiki/Measurable_function)
`f` that transforms values distributed according to measure `μ` to
values distributed according to a measure `ν`.

    y = vartransform(ν, μ, x)

Transforms a value `x` distributed according to `μ` to a value `y` distributed
according to `ν`.

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

* `MeasureBase.vartransform(ν::SomeStdMeasure, μ::CustomMeasure, x) = ...`
* `MeasureBase.vartransform(ν::MyMeasure, μ::SomeStdMeasure, x) = ...`

and/or

* `MeasureBase.vartransform_origin(ν::MyMeasure) = SomeMeasure(...)`
* `MeasureBase.from_origin(μ::MyMeasure, y) = x`
* `MeasureBase.to_origin(μ::MyMeasure, x) = y`

and ensure `MeasureBase.effndof(μ::MyMeasure)` is defined correctly.

If no direct transformation rule is available, `vartransform(ν, μ, x)` uses
the following strategy:

* Evaluate [`vartransform_origin`](@ref) for μ and ν. If both have an origin,
  select one as an intermediate measure using
  [`select_vartransform_intermediate`](@ref). Try to transform from `μ` to
  that intermediate measure and then to `ν` origin(s) of `μ` and/or `ν` if
  available.

* If all else fails, try to transform from μ to a standard multivariate
  uniform measure and then to ν.
"""
function vartransform end


function _vartransform_with_intermediate(ν, m, μ, x)
    x_m = vartransform(m, μ, x)
    _vartransform_with_intermediate_step2(ν, m, x_m)
end

@inline _vartransform_with_intermediate_step2(ν, m, x_m) = vartransform(ν, m, x_m)
@inline _vartransform_with_intermediate_step2(ν, m, x_m::NoTransformOrigin) = x_m

function _vartransform_with_intermediate(ν, m::NoTransformOrigin, μ, x)
    _vartransform_with_intermediate(ν, StdUniform()^effndof(μ), μ, x)
end


# Prevent endless recursion:
_vartransform_with_intermediate(::NU, ::NU, ::MU, x) where {NU,MU} = NoVarTransform{NU,MU}()
_vartransform_with_intermediate(::NU, ::MU, ::MU, x) where {NU,MU} = NoVarTransform{NU,MU}()

function vartransform(ν, μ, x)
    require_same_effndof(ν, μ)
    m = vartransform_intermediate(vartransform_origin(ν), vartransform_origin(μ))
    _vartransform_with_intermediate(ν, m, μ, x)
end

vartransform(::Any, ::Any, x::NoTransformOrigin) = x


"""
    struct VarTransformation <: Function

Transforms a variate from one measure to a variate of another.

In general users should not call `VarTransformation` directly, call
[`vartransform`](@ref) instead.
"""
struct VarTransformation{NU,MU} <: Function
    ν::NU
    μ::MU

    function VarTransformation{NU,MU}(ν::NU, μ::MU) where {NU,MU}
        require_same_effndof(ν, μ)
        return new{NU,MU}(ν, μ)
    end

    function VarTransformation(ν::NU, μ::MU) where {NU,MU}
        require_same_effndof(ν, μ)
        return new{NU,MU}(ν, μ)
    end
end

vartransform(ν, μ) = VarTransformation(ν, μ)


(f::VarTransformation)(x) = vartransform(f.ν, f.μ, x)

InverseFunctions.inverse(f::VarTransformation) = VarTransformation(f.μ, f.ν)


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





"""
    MeasureBase.select_vartransform_intermediate(a, b)

Selects one of two candidate pullback measures `a, b` to use as an
intermediate in variate transformations.

See [`MeasureBase.vartransform_intermediate`](@ref).
"""
function select_vartransform_intermediate end

select_vartransform_intermediate(nu, ::NoTransformOrigin) = nu
select_vartransform_intermediate(::NoTransformOrigin, mu) = mu
select_vartransform_intermediate(::NoTransformOrigin, mu::NoTransformOrigin) = mu

# Ensure forward and inverse transformation use the same intermediate:
@generated function select_vartransform_intermediate(a, b)
    return nameof(a) < nameof(b) ? :a : :b
end
