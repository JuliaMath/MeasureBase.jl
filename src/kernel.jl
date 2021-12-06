
export Kernel

abstract type AbstractKernel <: AbstractMeasure end

struct Kernel{F,G} <: AbstractKernel
    f::F
    g::G

    Kernel(::Type{T}, g::G) where {T,G} = new{Type{T},G}(T, g)
    Kernel(f::F, g::G) where {F,G} = new{F,G}(f, g)
end

struct ParameterizedKernel{F,N,T} <: AbstractKernel
    f::F
    param_maps::NamedTuple{N,T}

    ParameterizedKernel(::Type{F}, param_maps::NamedTuple{N,T}) where {F,N,T} =
        new{Type{F},N,T}(F, param_maps)
    ParameterizedKernel(f::F, param_maps::NamedTuple{N,T}) where {F,N,T} =
        new{F,N,T}(f, param_maps)
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

mapcall(t, x) = map(func -> func(x), t)

# (k::Kernel{Type{P},<:Tuple})(x) where {P<:ParameterizedMeasure} = k.f(mapcall(k.param_maps, x)...)

(k::ParameterizedKernel)(x) where {M} = k.f(; mapcall(k.param_maps, x)...)

(k::ParameterizedKernel)(x...) = k(x)

function (k::ParameterizedKernel)(x::Tuple)
    k.f(NamedTuple{k.param_maps}(x))
end

(k::Kernel)(x) = (k.f ∘ k.g)(x)

"""
For any `k::Kernel`, `basekernel` is expected to satisfy
```
basekernel(k)(p) == (basemeasure ∘ k)(p)
```

The main purpose of `basekernel` is to make it efficient to compute
```
basemeasure(d::ProductMeasure) = productmeasure(basekernel(d.f), d.xs)
```
"""
function basekernel end

# TODO: Find a way to do better than this
basekernel(f) = basemeasure ∘ f

basekernel(k::ParameterizedKernel) = kernel(basekernel(k.f), k.param_maps)

basekernel(k::Kernel) = kernel(basekernel(k.f), k.g)

basekernel(f::Returns) = Returns(basemeasure(f.value))

# export kernelize

# function kernelize(μ::M) where {N, M <: ParameterizedMeasure{N}}
#     C = constructorof(M)
#     (Kernel{C,}(NamedTuple{N}, ), values(getfield(μ, :par)))
# end

function Base.show(io::IO, μ::AbstractKernel)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

function Pretty.quoteof(k::Kernel)
    qf = Pretty.quoteof(k.f)
    qg = Pretty.quoteof(k.g)
    :(Kernel($qf, $qg))
end

function Pretty.quoteof(k::Kernel{F,typeof(identity)}) where {F}
    qf = Pretty.quoteof(k.f)
    :(Kernel($qf))
end

function Pretty.quoteof(k::ParameterizedKernel)
    qf = Pretty.quoteof(k.f)
    qg = Pretty.quoteof(k.param_maps)
    :(Kernel($qf, $qg))
end

export kernelfactor
