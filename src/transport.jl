"""
    MeasureBase.transport_origin(ν)

Default measure to pullback to resp. pushforward from when transforming
between `ν` and another measure.
"""
function transport_origin end

transport_origin(ν::NU) where {NU} = ν

"""
    MeasureBase.from_origin(ν, x)

Push `x` from `MeasureBase.transport_origin(μ)` forward to `ν`.
"""
function from_origin end

from_origin(ν::NU, x) where {NU} = x

"""
    MeasureBase.to_origin(ν, y)

Pull `y` from `ν` back to `MeasureBase.transport_origin(ν)`.
"""
function to_origin end

to_origin(ν::NU, x) where {NU} = x

"""
    struct MeasureBase.NoTransport{NU,MU} end

Indicates that no transformation from a measure of type `MU` to a measure of
type `NU` could be found.
"""
struct NoTransport{NU,MU,X} end

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

@inline transport_to(ν, μ) = TransportFunction(ν, μ)

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

function transport_def(m1::M, m2::M, x) where {M}
    @assert m1 === m2
    return x
end

export transport_rel

@inline function transport_rel(μ::M, ν::N, x::X) where {M,N,X}
    if static_hasmethod(transport_def, Tuple{M,N,X})
        return transport_def(μ, ν, x)
    end
    μs = origin_sequence(μ)
    νs = origin_sequence(ν)
    co = common_origin(μs, νs, X)
    isnothing(co) && begin
        μ = μs[end]
        ν = νs[end]
        @warn """
        No common transport origin for
            $μ
        and
            $ν

        To add one, add a method
            logdensity_def($μ, $ν, x)
        """
        return NoTransport{M,N,X}
    end
    return _transport_rel(μs, νs, co, x)
end


@generated function _transport_rel(
    νs::Tν,
    μs::Tμ,
    ::Tuple{StaticInt{N},StaticInt{M}},
    x::X,
) where {Tμ,Tν,M,N,X}
    sμ = schema(Tμ)
    sν = schema(Tν)
    
    xnames = map(m -> Symbol(:x_, m), 1:M)
    ynames = map(n -> Symbol(:y_, n), 1:N)

    q = quote
        $(Expr(:meta, :inline))
        x_1 = x
    end

    for j in 2:M
        x_old = xnames[j-1]
        x_new = xnames[j]
        push!(q.args, :($x_new = to_origin(μs[$(j-1)], $x_old)))
    end

    push!(q.args, :($(last(ynames)) = transport_def(last(νs), last(μs), $(last(xnames)))))

    for j in N-1:-1:1
        push!(q.args, :($(ynames[j]) = from_origin(νs[$j], $(ynames[j+1]))))
    end

    push!(q.args, :($(ynames[1])))
    return q
end



"""
    origin_sequence(m)

Construct the longest `Tuple` starting with `m` having each term as the base
measure of the previous term, and with no repeated entries.
"""
@inline function origin_sequence(μ::M) where {M}
    b_1 = μ
    done = false
    Base.Cartesian.@nexprs 10 i -> begin  # 10 is just some "big enough" number
        b_{i + 1} = if done
            nothing
        else
            transport_origin(b_{i})
        end
        if b_{i + 1} isa typeof(b_{i})
            done = true
            b_{i + 1} = nothing
        end
    end
    return filter(!isnothing, Base.Cartesian.@ntuple 10 b)
end


common_origin(μ, ν) = common_origin(μ, ν, Any)

"""
    common_origin(μ, ν, T) -> Tuple{StaticInt{i}, StaticInt{j}}

Find minimal (with respect to their sum) `i` and `j` such that there is a method

    logdensity_def(origin_sequence(μ)[i], origin_sequence(ν)[j], ::T)

This is used in `logdensity_rel` to help make that function efficient.
"""
@inline function common_origin(μ, ν, ::Type{T}) where {T}
    return common_origin(origin_sequence(μ), origin_sequence(ν), T)
end

@generated function common_origin(μ::M, ν::N, ::Type{T}) where {M<:Tuple,N<:Tuple,T}
    m = schema(M)
    n = schema(N)

    sols = Iterators.filter(
        ((i, j),) -> static_hasmethod(transport_def, Tuple{m[i],n[j],T}),
        Iterators.product(1:length(m), 1:length(n)),
    )
    isempty(sols) && return :(nothing)
    minsol = static.(argmin(((i, j),) -> i + j, sols))
    quote
        $minsol
    end
end
