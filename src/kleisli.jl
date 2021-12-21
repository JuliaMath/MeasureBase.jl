# TODO: Dangerous to export this - let's not
abstract type AbstractKleisli <: AbstractMeasure end

struct ParameterizedKleisli{F,N,T} <: AbstractKleisli
    f::F
    param_maps::NamedTuple{N,T}

    ParameterizedKleisli(::Type{F}, param_maps::NamedTuple{N,T}) where {F,N,T} =
        new{Type{F},N,T}(F, param_maps)
    ParameterizedKleisli(f::F, param_maps::NamedTuple{N,T}) where {F,N,T} =
        new{F,N,T}(f, param_maps)
end

"""
    kleisli(f, M)
    kleisli((f1, f2, ...), M)

A kleisli `κ = kleisli(f, m)` returns a wrapper around
a function `f` giving the parameters for a measure of type `M`,
such that `κ(x) = M(f(x)...)`
respective `κ(x) = M(f1(x), f2(x), ...)`

If the argument is a named tuple `(;a=f1, b=f1)`, `κ(x)` is defined as
`M(;a=f(x),b=g(x))`.

# Reference

* https://en.wikipedia.org/wiki/Markov_kleisli
"""
function kleisli end


# kleisli(Normal) do x
#     (μ=x,σ=x^2)
# end

kleisli(f, ::Type{M}) where {M} = kleisli(M, f)

mapcall(t, x) = map(func -> func(x), t)

# (k::Kleisli{Type{P},<:Tuple})(x) where {P<:ParameterizedMeasure} = k.f(mapcall(k.param_maps, x)...)

(k::ParameterizedKleisli)(x) where {M} = k.f(; mapcall(k.param_maps, x)...)

(k::ParameterizedKleisli)(x...) = k(x)

function (k::ParameterizedKleisli)(x::Tuple)
    k.f(NamedTuple{k.param_maps}(x))
end


"""
For any `k::Kleisli`, `basekleisli` is expected to satisfy
```
basekleisli(k)(p) == (basemeasure ∘ k)(p)
```

The main purpose of `basekleisli` is to make it efficient to compute
```
basemeasure(d::ProductMeasure) = productmeasure(basekleisli(d.f), d.xs)
```
"""
function basekleisli end

# TODO: Find a way to do better than this
basekleisli(f) = basemeasure ∘ f

basekleisli(k::ParameterizedKleisli) = kleisli(basekleisli(k.f), k.param_maps)

basekleisli(f::Returns) = Returns(basemeasure(f.value))


function Base.show(io::IO, μ::AbstractKleisli)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

function Pretty.quoteof(k::ParameterizedKleisli)
    qf = Pretty.quoteof(k.f)
    qg = Pretty.quoteof(k.param_maps)
    :(ParameterizedKleisli($qf, $qg))
end

