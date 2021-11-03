
export Kernel

abstract type AbstractKernel <: AbstractMeasure end

struct Kernel{F,S} <: AbstractKernel
    f::F
    ops::S

    Kernel(::Type{T}, ops::S) where {T,S} = new{Type{T},S}(T, ops)
    Kernel(f::F, ops::S) where {F,S} = new{F,S}(f, ops)
end

function Pretty.quoteof(k::Kernel)
    qf = Pretty.quoteof(k.f)
    qops = Pretty.quoteof(k.ops)
    :(Kernel($qf, $qops))
end

"""
    kernel(f, M)
    kernel((f1, f2, ...), M)

A kernel `κ = kernel(f, m)` returns a wrapper around
a function `f` giving the parameters for a measure of type `M`,
such that `κ(x) = M(f(x)...)`
respective `κ(x) = M(f1(x), f2(x), ...)`

If the argument is a named tuple `(;a=f1, b=f1)`, `κ(x)` is defined as
`M(;a=f(x),b=g(x))`.

# Reference

* https://en.wikipedia.org/wiki/Markov_kernel
"""
function kernel end

export kernel

# kernel(Normal) do x
#     (μ=x,σ=x^2)
# end

kernel(f, ::Type{M}) where {M} = kernel(M, f)

# TODO: Would this benefit from https://github.com/tisztamo/FunctionWranglers.jl?
mapcall(t, x) = map(func -> func(x), t)

# (k::Kernel{Type{P},<:Tuple})(x) where {P<:ParameterizedMeasure} = k.f(mapcall(k.ops, x)...)

(k::Kernel{M,<:NamedTuple})(x) where {M} = k.f(; mapcall(k.ops, x)...)

(k::Kernel{F,S})(x...) where {F,N,S<:NTuple{N,Symbol}} = k(x)

function (k::Kernel{F,S})(x::Tuple) where {F,N,S<:NTuple{N,Symbol}}
    k.f(NamedTuple{k.ops}(x))
end

(k::Kernel)(x) = k.f(k.ops(x))

"""
For any `k::Kernel`, `basekernel` is expected to satisfy
```
basemeasure(k(p)) == basekernel(k)(p)
```

The main purpose of `basekernel` is to make it efficient to compute
```
basemeasure(d::ProductMeasure) = productmeasure(basekernel(d.f), d.pars)
```
"""
function basekernel end

# TODO: Find a way to do better than this
basekernel(f) = basemeasure ∘ f

basekernel(k::Kernel) = kernel(basekernel(k.f), k.ops)
basekernel(f::Returns) = Returns(basemeasure(f.value))

# export kernelize

# function kernelize(μ::M) where {N, M <: ParameterizedMeasure{N}}
#     C = constructorof(M)
#     (Kernel{C,}(NamedTuple{N}, ), values(getfield(μ, :par)))
# end

export kernelfactor
