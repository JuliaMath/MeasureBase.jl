
"""
    abstract type ImplicitlyMapped

Supertype for objects that have been mapped in an implicit way.

The explicit map/function can only be determined given some kind of observed
result `obs` using

```julia
f_map = explicit_mapfunc(mapped::ImplicitlyMapped, obs)
```

The original object that has been implicitly mapped
may be retrieved via

```julia
obj = explicit_mapfunc(mapped::ImplicitlyMapped, obs)
```

Note that `obs` is typically *not* the directly result of `f_map(ob)`. Instead,
the relationship between `obj`, `f_map`, and `obs` depends on what `obj` is:

* A measure `mu = obj`: The mapping process is equivalent to
  `mapped_mu = pushfwd(f_map, mu, PushfwdRootMeasure())` and `obs` is an
  element of the measurable space of `mu`. Implicitly mapped measures support

  ```julia
  DensityInterface.DensityKind(mapped_mu::ImplicitlyMapped)
  DensityInterface.logdensityof(mapped_mu::ImplicitlyMapped, obs)
  ```

  and the explicitly mapped measure can be generated via

  ```julia
  explicit_measure(mapped_mu::ImplicitlyMapped, obs)
  ```

* A transition/Markov kernel `f_kernel = obj`, i.e. a function that maps
  points in some space to measures on a (possibly different) space:
  The mapping process is equivalent to
  `mapped_f_kernel = (p -> pushfwd(f_map, f_kernel(p), PushfwdRootMeasure()))`
  and `obs` is an element of the measurable space of the measures generated
  by the mapped kernel. Implicitly mapped transition/Markov kernels support

  ```julia
  Likelihood(mapped_f_kernel::ImplicitlyMapped, obs)
  ```

  and the explicitly mapped kernel can be generated via

  ```julia
  explicit_measure(mapped_mu::ImplicitlyMapped, obs)
  ```

# Implementation

Subtypes of `ImplicitlyMapped` that should support origin measures of
type `SomeRelevantMeasure` and observations of type `SomeRelevantObs`,
resulting in explicit maps/functions of type `SomeMapFunc`, must
implement/specialize

```julia
MeasureBase.implicit_origin(mapped::MyImplicitlyMapped)
MeasureBase.explicit_mapfunc(mapped::MyImplicitlyMapped, obs::SomeRelevantObs)::SomeMapFunc
```

and (except if functions of type `SomeMapFunc` are invertible via
`InverseFunctions.inverse`) must also specialize

```julia
MeasureBase.pushfwd(f::SomeMapFunc, mu::SomeRelevantMeasure, ::PushfwdRootMeasure)
```

Subtypes of `ImplicitlyMapped` may support multiple combinations of
observation and measure types.
"""
abstract type ImplicitlyMapped end
export ImplicitlyMapped

"""
    implicit_origin(mapped::ImplicitlyMapped)

Get the original object (a measure or transition/Markov kernel) that was
implicitly mapped.

See [ImplicitlyMapped](@ref) for detailed semantics.

# Implementation

`implicit_origin` must be implemented for subtypes of `ImplicitlyMapped`,
there is no default implementation.
"""
function implicit_origin end
export implicit_origin

"""
    explicit_mapfunc(mapped::ImplicitlyMapped, obs)

Get an explicit map/function based on an implicitly mapped object and an
observation.

See [ImplicitlyMapped](@ref) for detailed semantics.

# Implementation

`explicit_mapfunc` must be implemented for subtypes of `ImplicitlyMapped`,
there is no default implementation.
"""
function explicit_mapfunc end
export explicit_mapfunc

"""
    explicit_measure(mapped::ImplicitlyMapped, obs)

Get an explicitly mapped measure based on an implicitly mapped measure and an
observation that provides context on which pushforward to use on the onmapped
original measure `implicit_origin(mapped)`.

Used [`explicit_mapfunc`](@ref) to get the function to use in the pushforward.

# Implementation

`explicit_measure` does not need to be specialized for subtypes of
`ImplicitlyMapped`.
"""
function explicit_measure(mapped_measure::ImplicitlyMapped, obs)
    f_map = explicit_mapfunc(mapped_measure, obs)
    mu = implicit_origin(mapped_measure)
    return pushfwd(f_map, mu, PushfwdRootMeasure())
end
export explicit_measure

function DensityInterface.logdensityof(mapped_measure::ImplicitlyMapped, obs)
    return logdensityof(explicit_measure(mapped_measure, obs), obs)
end

function DensityInterface.DensityKind(mapped::ImplicitlyMapped)
    DensityKind(implicit_origin(mapped))
end

"""
    explicit_kernel(mapped::ImplicitlyMapped, obs)

Get an expliclity mapped transition/Markov kernel, based on an implicitly
mapped kernel and an observation that provides context on which pushforward
to add to the unmapped original kernel `implicit_origin(mapped)`.

Used [`explicit_mapfunc`](@ref) to get the function to use in the pushforward.

# Implementation

`explicit_kernel` does not need to be specialized for subtypes of
`ImplicitlyMapped`.
"""
function explicit_kernel(mapped_kernel::ImplicitlyMapped, obs)
    f_map = explicit_mapfunc(mapped_kernel, obs)
    f_kernel = implicit_origin(mapped_kernel)
    return (p -> pushfwd(f_map, f_kernel(p), PushfwdRootMeasure()))
end
export explicit_kernel

function Likelihood(mapped_kernel::ImplicitlyMapped, obs)
    return Likelihood(explicit_kernel(mapped_kernel, obs), obs)
end

"""
    struct MeasureBase.TakeAny{T} <: Function

Represents a function that takes n values from a collection.

`f = TakeAny(n)` treats all collections as unordered: `f(xs) may take the
first `n` elements of `xs`, but there is no guarantee. It must, however,
always take take the same elements from collections that are identical.

Constructor: `TakeAny(n::Union{Integer,Static.StaticInteger})`.
"""
struct TakeAny{T<:IntegerLike}
    n::T
end

_takeany_range(f::TakeAny, idxs) = first(idxs):first(idxs)+dynamic(f.n)-1
@inline _takeany_range(f::TakeAny, ::OneTo) = OneTo(dynamic(f.n))

@inline _takeany_range(::TakeAny{<:Static.StaticInteger{N}}, ::OneTo) where {N} = SOneTo(N)
@inline _takeany_range(::TakeAny{<:Static.StaticInteger{N}}, ::SOneTo) where {N} = SOneTo(N)

@inline (f::TakeAny)(xs::Tuple) = xs[begin:begin+f.n-1]
@inline (f::TakeAny)(xs::AbstractVector) = xs[_takeany_range(f, eachindex(xs))]

function (f::TakeAny)(xs)
    n = dynamic(f.n)
    ys = collect(Iterators.take(xs, n))
    length(ys) != n &&
        throw(ArgumentError("Can't take $n elements from a sequence shorter than $n"))
    return typeof(xs)(ys)
end

"""
    struct Marginalized{T} <: ImplicitlyMapped

Represents an implicitly marginalized measure or transition kernel.

Constructors:

* `Marginalized(mu)`
* `Marginalized(f_kernel)`

See [ImplicitlyMapped](@ref) for detailed semantics.

Example:

```julia
mu = productmeasure((a = StdUniform(), b = StdNormal(), c = StdExponential()))
obs = (a = 0.7, c = 1.2)

marg_mu_equiv = productmeasure((a = StdUniform(), c = StdExponential()))

logdensityof(Marginalized(mu), obs) ≈ logdensityof(marg_mu_equiv, obs)
```
"""
struct Marginalized{T} <: ImplicitlyMapped
    obj::T
end
export Marginalized

implicit_origin(mapped::Marginalized) = mapped.obj

function explicit_mapfunc(::Marginalized, obs::NamedTuple{names}) where {names}
    PropSelFunction{names,names}()
end
function pushfwd(f::PropSelFunction, mu::ProductMeasure{<:NamedTuple}, ::PushfwdRootMeasure)
    productmeasure(f(marginals(mu)))
end

explicit_mapfunc(::Marginalized, obs::AbstractVector) = TakeAny(length(obs))
explicit_mapfunc(::Marginalized, obs::StaticArray{Tuple{N}}) where {N} = TakeAny(static(N))

function pushfwd(
    f::TakeAny,
    mu::PowerMeasure{<:Any,<:Tuple{<:AbstractUnitRange}},
    ::PushfwdRootMeasure,
)
    n = f.n
    n_mu = length(mu)
    n_mu < n && throw(
        ArgumentError("Can't marginalize $n_mu dimensional power measure to $n dimensions"),
    )
    mu.parent^f.n
end
