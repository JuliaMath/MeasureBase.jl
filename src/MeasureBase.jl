module MeasureBase

const logtwo = log(2.0)

using Random
import Random: rand!

using FillArrays
using ConcreteStructs
using MLStyle

export ≪
export sampletype

export AbstractMeasure

abstract type AbstractMeasure end

import PrettyPrinting

const Pretty = PrettyPrinting

sampletype(μ::AbstractMeasure) = typeof(testvalue(μ))

# sampletype(μ::AbstractMeasure) = sampletype(basemeasure(μ))

export logdensity
export basemeasure
export basekernel

using LogExpFunctions: logsumexp

"""
    logdensity(μ::AbstractMeasure{X}, x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `basemeasure(μ)`.

Methods for computing density relative to other measures will be
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

include("kernel.jl")
include("parameterized.jl")
include("combinators/mapsto.jl")
include("combinators/half.jl")
include("exp.jl")
include("domains.jl")
include("utils.jl")
include("absolutecontinuity.jl")
include("macros.jl")

include("primitive.jl")
include("primitives/counting.jl")
include("primitives/lebesgue.jl")
include("primitives/dirac.jl")
include("primitives/trivial.jl")

include("densityinterface.jl")

include("combinators/factoredbase.jl")
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
