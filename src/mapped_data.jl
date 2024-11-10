"""
    MeasureBase.is_mapped_obs(::Type{Obs})

Return `Static.True()` or `Static.False()` if observations of type `Obs` may
be mapped observations (resp. a mapped points in a measurable space in
general).

A return value of `Static.True()` indicates that `map_result_and_func(obs::Obs)`
may return a different result than `(identity, obs)`. A return value of
`Static.False()` gurantees that `map_result_and_func(obs::Obs)` will always
return `(identity, obs)`.

See also [`map_result_and_func`](@ref).
"""
function is_mapped_obs end

# ToDo: Use generated function?
@inline function is_mapped_obs(::Type{T}) where {T}
    if T isa DataType
        for ST in T.parameters
            ST isa Type && is_mapped_obs(ST) isa True && return True()
        end
        return False()
    elseif T isa Union
        is_mapped_obs(T.a) isa True && return True()
        is_mapped_obs(T.b) isa True && return True()
        return False()
    elseif T isa UnionAll
        return is_mapped_obs(T.body)
    else
        # Shouldn't reach this point
        return False()
    end
end

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
f, y == map_result_and_func(WasMarginalized((a = 0.7, c = 1.2)))
x = (a = 0.7, b = something, c = 1.2, d = something_else)
````

will result in

```julia
y == (a = 0.7, c = 1.2)
y == f(x)
```

# Implementation

To when specializing `map_result_and_func` for additonal data types,
[`MeasureBase.is_mapped_obs`] must be specialized as well to return
`Static.True()` for the type.
"""
function map_result_and_func end
export map_result_and_func

#!!!!! ToDo: Imlement properly (needs a generated function)
@inline map_result_and_func(tagged_y) = (identity, tagged_y)

"""
    mapped_measure_and_obs(μ, maybe_tagged_y)

If `maybe_tagged_y` implies a mapping function `f_map` and a mapped value `y`,
construct the pushforward `ν` of `μ` under `f_map` and return `(f_map, y)`,
otherwise returns `(μ, tagged_y)`.

See [`map_result_and_func`](@ref) for details on mapped observations.
"""
function mapped_measure_and_obs(μ, tagged_y)
    _mapped_measure_and_obs_impl(is_mapped_obs(tagged_y), μ, tagged_y)
end
export mapped_measure_and_obs

function _mapped_measure_and_obs_impl(::True, μ, tagged_y)
    f_map, y = map_result_and_func(tagged_y)
    ν = pushfwd(f_map, μ)
    return (ν, y)
end

_mapped_measure_and_obs_impl(::False, μ, tagged_y) = (μ, tagged_y)

"""
    mapped_kernel_and_obs(f_kernel, maybe_tagged_y)

If `maybe_tagged_y` implies a mapping function `f_map` and a mapped value `y`,
construct `mapped_f_kernel` by composing a a pushforward under `f` with
`f_map` and returns `(mapped_f_kernel, y)`. Otherwise, return
`(f_kernel, maybe_tagged_y)`.

See [`map_result_and_func`](@ref) for details on mapped observations.
"""
@inline function mapped_kernel_and_obs(f_kernel, tagged_y)
    _mapped_kernel_and_obs_impl(is_mapped_obs(tagged_y), f_kernel, tagged_y)
end
export mapped_kernel_and_obs

function _mapped_kernel_and_obs_impl(::True, f_kernel, tagged_y)
    f_map, y = map_result_and_func(tagged_y)
    mapped_k = _map_kernel(f_map, f_kernel)

    return (mapped_k, y)
end
_map_kernel(f_map, f_kernel) = pushfwd(f_map, PushfwdRootMeasure()) ∘ f_kernel
_map_kernel(::typeof(identity), f_kernel) = f_kernel

_mapped_kernel_and_obs_impl(::False, f_kernel, tagged_y) = (f_kernel, tagged_y)

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
    struct WasMarginalized{T}

Annotates a data/observation object to indicate that is was the result of
marginalization.

Constructors:

* `WasMarginalized(x)`

See [`map_result_and_func`](@ref) for more details.
"""
struct WasMarginalized{T} <: MappedObservation
    obj::T
end
export WasMarginalized

implicit_origin(mapped::WasMarginalized) = mapped.obj

function explicit_mapfunc(::WasMarginalized, obs::NamedTuple{names}) where {names}
    PropSelFunction{names,names}()
end
function pushfwd(f::PropSelFunction, mu::ProductMeasure{<:NamedTuple}, ::PushfwdRootMeasure)
    productmeasure(f(marginals(mu)))
end

explicit_mapfunc(::WasMarginalized, obs::AbstractVector) = TakeAny(length(obs))
function explicit_mapfunc(::WasMarginalized, obs::StaticArray{Tuple{N}}) where {N}
    TakeAny(static(N))
end

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
