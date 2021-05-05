module MeasureBase


using Random

using ConcreteStructs
using MLStyle
using KeywordCalls
using Compat

export ≪
export sampletype

export AbstractMeasure

abstract type AbstractMeasure end

sampletype(μ::AbstractMeasure) = typeof(testvalue(μ))

# sampletype(μ::AbstractMeasure) = sampletype(basemeasure(μ))

export logdensity

"""
    logdensity(μ::AbstractMeasure{X}, x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `basemeasure(μ)`.

Methods for computing density relative to other measures will be
"""
function logdensity end

include("exp.jl")
include("domains.jl")
include("utils.jl")
include("absolutecontinuity.jl")
include("basemeasures.jl")
include("parameterized.jl")
include("macros.jl")
include("combinators/weighted.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("combinators/for.jl")
include("combinators/power.jl")
include("combinators/likelihood.jl")
include("combinators/elementwise.jl")
include("combinators/spikemixture.jl")
include("rand.jl")
include("density.jl")
# include("pushforward.jl")
include("kernel.jl")

end
