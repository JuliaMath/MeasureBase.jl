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

import IfElse: ifelse
export logdensity_def
export basemeasure
export basekleisli

"""
    inssupport(m, x)
    insupport(m)

`insupport(m,x)` computes whether `x` is in the support of `m`.

`insupport(m)` returns a function, and satisfies

    insupport(m)(x) == insupport(m, x)
"""
function insupport end

export insupport

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

using LogExpFunctions: logsumexp

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

include("schema.jl")
include("splat.jl")
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

include("combinators/conditional.jl")
include("combinators/bind.jl")
include("combinators/transformedmeasure.jl")
include("combinators/factoredbase.jl")
include("combinators/weighted.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("combinators/power.jl")
include("combinators/spikemixture.jl")
include("combinators/likelihood.jl")
include("combinators/pointwise.jl")
include("combinators/restricted.jl")
include("combinators/smart-constructors.jl")

include("rand.jl")

include("density.jl")

include("interface.jl")

using .Interface

end # module MeasureBase