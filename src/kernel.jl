export AbstractTransitionKernel,
    GenericTransitionKernel, TypedTransitionKernel, ParameterizedTransitionKernel

abstract type AbstractTransitionKernel <: AbstractMeasure end

struct GenericTransitionKernel{F} <: AbstractTransitionKernel
    f::F
end

(k::GenericTransitionKernel)(x) = k.f(x)

struct TypedTransitionKernel{M,F} <: AbstractTransitionKernel
    m::M
    f::F
end

(k::TypedTransitionKernel)(x) = (k.m ∘ k.f)(x)
struct ParameterizedTransitionKernel{M,S,N,T} <: AbstractTransitionKernel
    m::M
    suff::S
    param_maps::NamedTuple{N,T}

    function ParameterizedTransitionKernel(
        ::Type{M},
        suff::S,
        param_maps::NamedTuple{N,T},
    ) where {M,S,N,T}
        new{Type{M},S,N,T}(M, suff, param_maps)
    end
    function ParameterizedTransitionKernel(
        m::M,
        suff::S,
        param_maps::NamedTuple{N,T},
    ) where {M,S,N,T}
        new{M,S,N,T}(m, suff, param_maps)
    end
end

"""
A *kernel* is a function that returns a measure.

    k1 = kernel() do x
        Normal(x, x^2)
    end

    k2 = kernel(Normal) do x
        (μ = x, σ = x^2)
    end

    k3 = kernel(Normal; μ = identity, σ = abs2)

    k4 = kernel(Normal; μ = first, σ = last) do x
        (x, x^2)
    end

    x = randn(); k1(x) == k2(x) == k3(x) == k4(x)

This function is not exported, because "kernel" can have so many other meanings.
See for example https://github.com/JuliaGaussianProcesses/KernelFunctions.jl for
another common use of this term.

# Reference

* https://en.wikipedia.org/wiki/Markov_kernel
"""
function kernel end

mapcall(t, x) = map(func -> func(x), t)

function (k::ParameterizedTransitionKernel)(x)
    s = k.suff(x)
    k.m(; mapcall(k.param_maps, s)...)
end

(k::AbstractTransitionKernel)(x1, x2, xs...) = k((x1, x2, xs...))

(k::AbstractTransitionKernel)(; kwargs...) = k(NamedTuple(kwargs))

"""
For any `k::TransitionKernel`, `basekernel` is expected to satisfy
```
basekernel(k)(p) == (basemeasure ∘ k)(p)
```

The main purpose of `basekernel` is to make it efficient to compute
```
basemeasure(d::ProductMeasure) == productmeasure(basekernel(d.f), d.xs)
```
"""
function basekernel end

# TODO: Find a way to do better than this
basekernel(f) = basemeasure ∘ f

basekernel(f::Returns) = Returns(basemeasure(f.value))

function Base.show(io::IO, μ::AbstractTransitionKernel)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

function Pretty.tile(k::K) where {K<:AbstractTransitionKernel}
    Pretty.list_layout(
        Pretty.tile.([getproperty(k, p) for p in propertynames(k)]),
        prefix = nameof(constructorof(K)),
    )
end

const kleisli = kernel

export kleisli

kernel(k::AbstractTransitionKernel) = k
