# TODO: Dangerous to export this - let's not
abstract type AbstractTransitionKernel <: AbstractMeasure end

struct ParameterizedTransitionKernel{F,N,T} <: AbstractTransitionKernel
    f::F
    param_maps::NamedTuple{N,T}

    ParameterizedTransitionKernel(::Type{F}, param_maps::NamedTuple{N,T}) where {F,N,T} =
        new{Type{F},N,T}(F, param_maps)
    ParameterizedTransitionKernel(f::F, param_maps::NamedTuple{N,T}) where {F,N,T} =
        new{F,N,T}(f, param_maps)
end

"""
    kernel(f, M)
    kernel((f1, f2, ...), M)

A kernel `κ = kernel(f, m)` returns a wrapper around a function `f` giving the
parameters for a measure of type `M`, such that `κ(x) = M(f(x)...)` respective
`κ(x) = M(f1(x), f2(x), ...)`

If the argument is a named tuple `(;a=f1, b=f1)`, `κ(x)` is defined as
`M(;a=f(x),b=g(x))`.

This function is not exported, because "kernel" can have so many other meanings.
See for example https://github.com/JuliaGaussianProcesses/KernelFunctions.jl for
another common use of this term.

# Reference

* https://en.wikipedia.org/wiki/Markov_kernel
"""
function kernel end


# kernel(Normal) do x
#     (μ=x,σ=x^2)
# end

kernel(f, ::Type{M}) where {M} = kernel(M, f)

mapcall(t, x) = map(func -> func(x), t)

# (k::TransitionKernel{Type{P},<:Tuple})(x) where {P<:ParameterizedMeasure} = k.f(mapcall(k.param_maps, x)...)

(k::ParameterizedTransitionKernel)(x) = k.f(; mapcall(k.param_maps, x)...)

(k::ParameterizedTransitionKernel)(x...) = k(x)

function (k::ParameterizedTransitionKernel)(x::Tuple)
    k.f(NamedTuple{k.param_maps}(x))
end


"""
For any `k::TransitionKernel`, `basekernel` is expected to satisfy
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

basekernel(k::ParameterizedTransitionKernel) = kernel(basekernel(k.f), k.param_maps)

basekernel(f::Returns) = Returns(basemeasure(f.value))


function Base.show(io::IO, μ::AbstractTransitionKernel)
    io = IOContext(io, :compact => true)
    Pretty.pprint(io, μ)
end

function Pretty.quoteof(k::ParameterizedTransitionKernel)
    qf = Pretty.quoteof(k.f)
    qg = Pretty.quoteof(k.param_maps)
    :(ParameterizedTransitionKernel($qf, $qg))
end

const kleisli = kernel

export kleisli