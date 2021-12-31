module MeasureBase

const logtwo = log(2.0)

using Random
import Random: rand!
import Random: gentype
using Statistics
using LinearAlgebra

import DensityInterface: logdensityof
import DensityInterface: densityof
import DensityInterface: DensityKind
using DensityInterface

import ConstructionBase
using ConstructionBase: constructorof

using PrettyPrinting
const Pretty = PrettyPrinting

using FillArrays
using Static

export ≪
export gentype
export rebase

export AbstractMeasure

abstract type AbstractMeasure end

using Static: @constprop

function Pretty.tile(d::M) where {M<:AbstractMeasure}
    the_names = fieldnames(typeof(d))
    result = Pretty.literal(repr(M))
    isempty(the_names) && return result * Pretty.literal("()")
    Pretty.list_layout(Pretty.tile.([getfield(d, n) for n in the_names]); prefix=result)
end

@inline DensityKind(::AbstractMeasure) = HasDensity()

gentype(μ::AbstractMeasure) = typeof(testvalue(μ))

# gentype(μ::AbstractMeasure) = gentype(basemeasure(μ))

export logdensity_def
export basemeasure
export basekleisli

using LogExpFunctions: logsumexp

"""
    logdensity_def(μ::AbstractMeasure{X}, x::X)

Compute the logdensity of the measure μ at the point x. This is the standard way
to define `logdensity` for a new measure. the base measure is implicit here, and
is understood to be `basemeasure(μ)`.

Methods for computing density relative to other measures will be
"""
function logdensity_def end

using Compat

include("proxies.jl")
include("kleisli.jl")
include("parameterized.jl")
include("combinators/half.jl")
include("domains.jl")
include("primitive.jl")
include("utils.jl")
include("absolutecontinuity.jl")

include("primitives/counting.jl")
include("primitives/lebesgue.jl")
include("primitives/dirac.jl")
include("primitives/trivial.jl")

include("combinators/bind.jl")
include("combinators/transformedmeasure.jl")
include("combinators/factoredbase.jl")
include("combinators/weighted.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("combinators/for.jl")
include("combinators/power.jl")
include("combinators/spikemixture.jl")
include("combinators/likelihood.jl")
include("combinators/pointwise.jl")
include("combinators/restricted.jl")
include("combinators/smart-constructors.jl")
include("combinators/conditional.jl")

include("rand.jl")

include("density.jl")
module Interface

using Reexport
using MeasureBase
using MeasureBase:basemeasure_depth, proxy
@reexport using Test

export test_interface
export basemeasure_depth
export proxy

include("help.jl")
include("interface.jl")
end # module Interface

using .Interface

end # module MeasureBase