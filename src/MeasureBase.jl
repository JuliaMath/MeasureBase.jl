module MeasureBase

using Base: @propagate_inbounds

using Random
import Random: rand!
import Random: gentype
using Statistics
using LinearAlgebra

import IntervalSets
# This seems harder than it should be to get `IntervalSets.:(..)`
@eval (using IntervalSets: $(Symbol(IntervalSets.:(..))))

using IntervalSets: Interval, width

import DensityInterface: logdensityof
import DensityInterface: densityof
import DensityInterface: DensityKind
using DensityInterface

using InverseFunctions
using ChangesOfVariables

import Base.iterate
import ConstructionBase
using ConstructionBase: constructorof

using PrettyPrinting
const Pretty = PrettyPrinting

using ChainRulesCore
using FillArrays
using Static

export ≪
export gentype
export rebase

export AbstractMeasure

import IfElse: ifelse
export logdensity_def
export basemeasure
export basekernel
export productmeasure

"""
    inssupport(m, x)
    insupport(m)

`insupport(m,x)` computes whether `x` is in the support of `m`.

`insupport(m)` returns a function, and satisfies

    insupport(m)(x) == insupport(m, x)
"""
function insupport end

export insupport
export getdof
export transport_to

include("insupport.jl")

abstract type AbstractMeasure end

AbstractMeasure(m::AbstractMeasure) = m

using Static: @constprop

function Pretty.quoteof(d::M) where {M<:AbstractMeasure}
    the_names = fieldnames(typeof(d))
    :($M($([getfield(d, n) for n in the_names]...)))
end

@inline DensityKind(::AbstractMeasure) = HasDensity()

Broadcast.broadcastable(m::AbstractMeasure) = Ref(m)

gentype(μ::AbstractMeasure) = typeof(testvalue(μ))

# gentype(μ::AbstractMeasure) = gentype(basemeasure(μ))

using NaNMath
using LogExpFunctions: logsumexp, logistic, logit

@deprecate instance_type(x) Core.Typeof(x) false

"""
`logdensity_def` is the standard way to define a log-density for a new measure.
Note that this definition does not include checking for membership in the
support; this is instead checked using `insupport`. `logdensity_def` is
a low-level function, and should typically not be called directly. See
`logdensityof` for more information and other alternatives.

---

    logdensity_def(m, x)

Compute the log-density of the measure m at the point `x`, relative to
`basemeasure(m)`, and assuming `insupport(m, x)`.

---

    logdensity_def(m1, m2, x)

Compute the log-density of `m1` relative to `m2` at the point `x`, assuming
`insupport(m1, x)` and `insupport(m2, x)`.
"""
function logdensity_def end

using Compat

using IrrationalConstants

include("getdof.jl")
include("transport.jl")
include("schema.jl")
include("splat.jl")
include("proxies.jl")
include("kernel.jl")
include("parameterized.jl")
include("combinators/half.jl")
include("domains.jl")
include("primitive.jl")
include("utils.jl")
include("mass-interface.jl")
# include("absolutecontinuity.jl")

include("primitives/counting.jl")
include("primitives/lebesgue.jl")
include("primitives/dirac.jl")
include("primitives/trivial.jl")

include("combinators/bind.jl")
include("combinators/transformedmeasure.jl")
include("combinators/weighted.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("combinators/power.jl")
include("combinators/spikemixture.jl")
include("combinators/likelihood.jl")
include("combinators/pointwise.jl")
include("combinators/restricted.jl")
include("combinators/smart-constructors.jl")
include("combinators/powerweighted.jl")
include("combinators/conditional.jl")

include("standard/stdmeasure.jl")
include("standard/stduniform.jl")
include("standard/stdexponential.jl")
include("standard/stdlogistic.jl")
include("standard/stdnormal.jl")

include("rand.jl")
include("fixedrng.jl")

include("density.jl")
include("density-core.jl")

include("interface.jl")

using .Interface

end # module MeasureBase
