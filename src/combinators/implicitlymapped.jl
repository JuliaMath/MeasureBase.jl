
"""
    abstract type MappedObservation

Supertype for observation/data that has been mapped, and so implies a
pushforward on measures in the context of `logdensityof` and `likelihoodof`.

The explicit map/function can only be determined given some kind of observed
result `obs` using

```julia
f_map = explicit_mapfunc(mapped::MappedObservation, obs)
```

The original object that has been implicitly mapped
may be retrieved via

```julia
obj = explicit_mapfunc(mapped::MappedObservation, obs)
```

Note that `obs` is typically *not* the directly result of `f_map(ob)`. Instead,
the relationship between `obj`, `f_map`, and `obs` depends on what `obj` is:

* A measure `mu = obj`: The mapping process is equivalent to
  `mapped_mu = pushfwd(f_map, mu, PushfwdRootMeasure())` and `obs` is an
  element of the measurable space of `mu`. Implicitly mapped measures support

  ```julia
  DensityInterface.DensityKind(mapped_mu::MappedObservation)
  DensityInterface.logdensityof(mapped_mu::MappedObservation, obs)
  ```

  and the explicitly mapped measure can be generated via

  ```julia
  explicit_measure(mapped_mu::MappedObservation, obs)
  ```

* A transition/Markov kernel `f_kernel = obj`, i.e. a function that maps
  points in some space to measures on a (possibly different) space:
  The mapping process is equivalent to
  `mapped_f_kernel = (p -> pushfwd(f_map, f_kernel(p), PushfwdRootMeasure()))`
  and `obs` is an element of the measurable space of the measures generated
  by the mapped kernel. Implicitly mapped transition/Markov kernels support

  ```julia
  Likelihood(mapped_f_kernel::MappedObservation, obs)
  ```

  and the explicitly mapped kernel can be generated via

  ```julia
  explicit_measure(mapped_mu::MappedObservation, obs)
  ```

# Implementation

Subtypes of `MappedObservation` that should support origin measures of
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

Subtypes of `MappedObservation` may support multiple combinations of
observation and measure types.
"""
abstract type MappedObservation end
export MappedObservation

"""
    map_result_and_func(tagged_y)

Get the map function associated with the a mapped observation and the
result of the maps.

Assuming `y` the result of a function `f(y)` and `tagged_y` provides enough
information to determine `f` any `y`, then

```
(f, y) == map_result_and_func(tagged_y)
```

`tagged_y` may be identical to y - if so, then `f` must be idempotent.

For example, with
mu = productmeasure((a = StdUniform(), b = StdNormal(), c = StdExponential()))

```julia
f, y == map_result_and_func(was_marginalized((a = 0.7, c = 1.2)))
x = (a = 0.7, b = something, c = 1.2, d = something_else)
````

will result in

```julia
y == (a = 0.7, c = 1.2)
y == f(x)
```
"""
function map_result_and_func end
export map_result_and_func

function DensityInterface.logdensityof(f_kernel::Base.Callable, x)
    f_map, x_unwrapped = map_result_and_func(x)
    mapped_kernel = _map_kernel(f_map, f_kernel)

    return (mapped_kernel, x_unwrapped)
end

function likelihoodof(f_kernel::Base.Callable, x)
    f_map, x_unwrapped = map_result_and_func(x)
    mapped_kernel = _map_kernel(f_map, f_kernel)

    return Likelihood{Base.Typeof(mapped_kernel),Base.Typeof(x_unwrapped)}(
        mapped_kernel,
        x_unwrapped,
    )
end

#!!!!!!!! make an object for this instead of using an anonymous function:
_map_kernel(f_map, f_kernel) = (p -> pushfwd(f_map, f_kernel(p), PushfwdRootMeasure()))
_map_kernel(::typeof(identity), f_kernel) = f_kernel

function DensityInterface.DensityKind(mapped::MappedObservation)
    DensityKind(implicit_origin(mapped))
end

"""
    explicit_kernel(mapped::MappedObservation, obs)

Get an expliclity mapped transition/Markov kernel, based on an implicitly
mapped kernel and an observation that provides context on which pushforward
to add to the unmapped original kernel `implicit_origin(mapped)`.

Used [`explicit_mapfunc`](@ref) to get the function to use in the pushforward.

# Implementation

`explicit_kernel` does not need to be specialized for subtypes of
`MappedObservation`.
"""
function explicit_kernel(mapped_kernel::MappedObservation, obs)
    f_map = explicit_mapfunc(mapped_kernel, obs)
    f_kernel = implicit_origin(mapped_kernel)
    return (p -> pushfwd(f_map, f_kernel(p), PushfwdRootMeasure()))
end
export explicit_kernel

function Likelihood(mapped_kernel::MappedObservation, obs)
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
    struct Marginalized{T} <: MappedObservation

Represents an implicitly marginalized measure or transition kernel.

Constructors:

* `Marginalized(mu)`
* `Marginalized(f_kernel)`

See [MappedObservation](@ref) for detailed semantics.

Example:

```julia
mu = productmeasure((a = StdUniform(), b = StdNormal(), c = StdExponential()))
obs = (a = 0.7, c = 1.2)

marg_mu_equiv = productmeasure((a = StdUniform(), c = StdExponential()))

logdensityof(Marginalized(mu), obs) ≈ logdensityof(marg_mu_equiv, obs)
```
"""
struct Marginalized{T} <: MappedObservation
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
