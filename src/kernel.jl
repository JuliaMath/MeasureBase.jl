"""Abstract supertype for all kernels."""
abstract type AbstractKernel <: AbstractMeasure end

"""
    Kernel{F,S} <: AbstractKernel

Store a kernel, i.e. a function `k` of an argument `x` that returns a measure `k(x)`.

# Fields
- `f::F`: type (& constructor) for the generated measure
- `ops::S`: typically a `NamedTuple` of operations that convert `x` into arguments for the measure constructor `f`

# Example

To define a Normal approximation to the Poisson distribution with parameter `x`, one would write
```julia
k = x -> Normal(μ=x, σ²=x)
```
The `Kernel` constructor stores this in a way that is easier for the Julia compiler to make sense of:
```julia
k = Kernel(Normal, (μ=identity, σ²=identity))
````

# Reference

- https://en.wikipedia.org/wiki/Markov_kernel
"""
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

# TODO: Would this benefit from https://github.com/tisztamo/FunctionWranglers.jl?
"""
    mapcall(t, x)

Apply each function in `t` to each component of `x`.

# Arguments
- `t`: tuple of functions
- `x`: argument(s) for the function(s) in `t`
"""
mapcall(t, x) = map(func -> func(x), t)

# (k::Kernel{Type{P},<:Tuple})(x) where {P<:ParameterizedMeasure} = k.f(mapcall(k.ops, x)...)

(k::Kernel{M,<:NamedTuple})(x) where {M} = k.f(; mapcall(k.ops, x)...)

(k::Kernel{F,S})(x...) where {F,N,S<:NTuple{N,Symbol}} = k(x)

function (k::Kernel{F,S})(x::Tuple) where {F,N,S<:NTuple{N,Symbol}}
    k.f(NamedTuple{k.ops}(x))
end

(k::Kernel)(x) = k.f(k.ops(x))

"""
    kernel(f, M)
    kernel((f1, f2, ...), M)

Create a [`Kernel`](@ref) `k` around a function `f` that gives the parameters for a measure of type `M`.

When `k` is called with an argument `x`, it will create a measure `k(x) = M(f(x)...)`, resp. `k(x) = M(f1(x), f2(x), ...)`
If the argument to `kernel` is a named tuple of functions `(a=f1, b=f2)`, then `k(x)` is defined as `M(;a=f1(x),b=f2(x))`.
"""
function kernel end

# kernel(Normal) do x
#     (μ=x,σ=x^2)
# end

kernel(f, ::Type{M}) where {M} = kernel(M, f)

"""
    basekernel(k::Kernel)

For any `k::Kernel`, `basekernel` is expected to satisfy
```
basemeasure(k(x)) == basekernel(k)(x)
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

# function kernelize(μ::M) where {N, M <: ParameterizedMeasure{N}}
#     C = constructorof(M)
#     (Kernel{C,}(NamedTuple{N}, ), values(getfield(μ, :par)))
# end
