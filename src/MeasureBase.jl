module MeasureBase

const logtwo = log(2.0)

using Random

using ConcreteStructs
using MLStyle

export ≪
export sampletype

export AbstractMeasure

abstract type AbstractMeasure end

sampletype(μ::AbstractMeasure) = typeof(testvalue(μ))

# sampletype(μ::AbstractMeasure) = sampletype(basemeasure(μ))

export logdensity
export basemeasure

using LogExpFunctions: logsumexp

"""
    logdensity(μ::AbstractMeasure{X}, x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `basemeasure(μ)`.

Methods for computing density relative to other measures will be
"""
function logdensity end

include("combinators/half.jl")
include("exp.jl")
include("domains.jl")
include("utils.jl")
include("absolutecontinuity.jl")
include("parameterized.jl")
include("macros.jl")
include("resettablerng.jl")

include("primitive.jl")
include("primitives/counting.jl")
include("primitives/lebesgue.jl")
include("primitives/dirac.jl")
include("primitives/trivial.jl")

include("combinators/factoredbase.jl")
include("combinators/weighted.jl")
include("combinators/affine.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("combinators/for.jl")
include("combinators/power.jl")
include("combinators/spikemixture.jl")
include("kernel.jl")
include("combinators/likelihood.jl")
include("combinators/pointwise.jl")

include("rand.jl")

include("density.jl")

end
