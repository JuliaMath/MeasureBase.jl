# TODO: Dangerous to export this - let's not
abstract type AbstractTransitionKernel <: AbstractMeasure end

struct GenericTransitionKernel{F} <: AbstractTransitionKernel
    f::F
end

(k::GenericTransitionKernel)(x) = k.f(x)
(k::GenericTransitionKernel)(x1, x2, xs...) = k.f((x1, x2, xs...))

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

# """
#     kernel(f, M)
#     kernel((f1, f2, ...), M)

# A kernel `κ = kernel(f, m)` returns a wrapper around a function `f` giving the
# parameters for a measure of type `M`, such that `κ(x) = M(f(x)...)` respective
# `κ(x) = M(f1(x), f2(x), ...)`

# If the argument is a named tuple `(;a=f1, b=f1)`, `κ(x)` is defined as
# `M(;a=f(x),b=g(x))`.

# This function is not exported, because "kernel" can have so many other meanings.
# See for example https://github.com/JuliaGaussianProcesses/KernelFunctions.jl for
# another common use of this term.

# # Reference

# * https://en.wikipedia.org/wiki/Markov_kernel
# """
# function kernel end

# # kernel(Normal) do x
# #     (μ=x,σ=x^2)
# # end

# function _kernel(::Type{M}, f::F, ::Type{T}) where {M,F<:Function,T}
#     TypedTransitionKernel{Type{M},F}(M,f)
# end

mapcall(t, x) = map(func -> func(x), t)

# # (k::TransitionKernel{Type{P},<:Tuple})(x) where {P<:ParameterizedMeasure} = k.f(mapcall(k.param_maps, x)...)

function (k::ParameterizedTransitionKernel)(x)
    s = k.suff(x)
    k.m(; mapcall(k.param_maps, s)...)
end

(k::AbstractTransitionKernel)(x1, x2, xs...) = k((x1, x2, xs...))

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

basekernel(f::Returns) = Returns(basemeasure(f.value))

# function Base.show(io::IO, μ::AbstractTransitionKernel)
#     io = IOContext(io, :compact => true)
#     Pretty.pprint(io, μ)
# end

# function Pretty.quoteof(k::ParameterizedTransitionKernel)
#     qf = Pretty.quoteof(k.f)
#     qg = Pretty.quoteof(k.param_maps)
#     :(ParameterizedTransitionKernel($qf, $qg))
# end

const kleisli = kernel

export kleisli

# kernel(f, pars::NamedTuple) = ParameterizedTransitionKernel(f, pars)

# # kernel(Normal{(:μ,), Tuple{Int64}})
# function kernel(::Type{M}) where {M<:AbstractMeasure}
#     constructorof(M)
# end

# # kernel(::Type{P}, op::O) where {O, N, P<:ParameterizedMeasure{N}} = kernel{constructorof(P),O}(op)

# function kernel(::Type{M}; param_maps...) where {M}
#     nt = NamedTuple(param_maps)
#     kernel(M, nt)
# end


kernel(k::AbstractTransitionKernel) = k