module MeasureBase

const logtwo = log(2.0)

using Random
import Random: rand!
import Random: gentype
import DensityInterface: logdensityof
import DensityInterface: densityof
import DensityInterface: DensityKind
using DensityInterface

using PrettyPrinting
const Pretty = PrettyPrinting

using FillArrays
using Static

export ≪
export gentype

export AbstractMeasure

abstract type AbstractMeasure end

@inline DensityKind(::AbstractMeasure) = HasDensity()

gentype(μ::AbstractMeasure) = typeof(testvalue(μ))


# gentype(μ::AbstractMeasure) = gentype(basemeasure(μ))

export logdensity_def
export basemeasure
export basekernel

using LogExpFunctions: logsumexp

"""
    logdensity_def(μ::AbstractMeasure{X}, x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `basemeasure(μ)`.

Methods for computing density relative to other measures will be
"""
function logdensity_def end

if VERSION < v"1.7.0-beta1.0"
    @eval begin
        struct Returns{T}
            value::T
        end

        (f::Returns)(x) = f.value
    end
end

include("proxies.jl")
include("kernel.jl")
include("parameterized.jl")
include("combinators/mapsto.jl")
include("combinators/half.jl")
include("domains.jl")
include("primitive.jl")
include("utils.jl")
include("absolutecontinuity.jl")

include("primitives/counting.jl")
include("primitives/lebesgue.jl")
include("primitives/dirac.jl")
include("primitives/trivial.jl")

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
