module MeasureBase

## Imports

using ConcreteStructs
using ConstructionBase: ConstructionBase
using FillArrays
using KeywordCalls
using LogExpFunctions: logsumexp
using MLStyle
using PrettyPrinting: PrettyPrinting
using Random: Random, AbstractRNG, rand!
using Tricks

## Exports

export AbstractMeasure
export Kernel

export ≪
export basemeasure
export sampletype
export logdensity

export basekernel
export kernel
export kernelfactor
# export kernelize

export ParameterizedMeasure
export params
export paramnames

export ℝ
export ℝ₊
export 𝕀
export ℤ
export IntegerRange

export testvalue
export rootmeasure

export @parameterized, @half

export ↦, mapsto

export FactoredBase

export Half

## Constants

const logtwo = log(2.0)

const Pretty = PrettyPrinting

## Basic types and functions

"""The abstract supertype of all measures."""
abstract type AbstractMeasure end

"""
    sampletype(μ::AbstractMeasure)

Return the type of a sample drawn from measure `μ`.
"""
sampletype(μ::AbstractMeasure) = typeof(testvalue(μ))

# sampletype(μ::AbstractMeasure) = sampletype(basemeasure(μ))

"""
    logdensity(μ::AbstractMeasure{X}, x::X)
    logdensity(μ, ν, x)

Compute the logdensity of the measure `μ` at the point `x`.

If `ν` is specified, the density is computed wrt `ν`. Otherwise, it is computed wrt the (local) base measure `basemeasure(μ, x)`.
This is the standard way to define `logdensity` for a new measure.
"""
function logdensity end

if VERSION < v"1.7.0-beta1.0"
    @eval begin
        struct Returns{T}
            value::T
        end

        (f::Returns)(x) = f.value
    end
end

## Includes

include("kernel.jl")
include("parameterized.jl")

include("exp.jl")
include("domains.jl")
include("utils.jl")
# include("absolutecontinuity.jl")
include("macros.jl")

include("primitive.jl")
include("primitives/counting.jl")
include("primitives/lebesgue.jl")
include("primitives/dirac.jl")
include("primitives/trivial.jl")

include("combinators/mapsto.jl")
include("combinators/factoredbase.jl")
include("combinators/half.jl")
include("combinators/weighted.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("combinators/for.jl")
include("combinators/power.jl")
include("combinators/affine.jl")
include("combinators/spikemixture.jl")
include("combinators/likelihood.jl")
include("combinators/pointwise.jl")
include("combinators/restricted.jl")
include("combinators/smart-constructors.jl")

include("rand.jl")

include("density.jl")

end
