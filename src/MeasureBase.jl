module MeasureBase

using Base: @propagate_inbounds
using Base: OneTo

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
using DensityInterface: FuncDensity, LogFuncDensity
using DensityInterface

using InverseFunctions
using InverseFunctions: FunctionWithInverse
using ChangesOfVariables
using ConstantRNGs

import Base.iterate
import ConstructionBase
using ConstructionBase: constructorof
using IntervalSets

import StaticArrays
using StaticArrays:
    StaticArray, StaticVector, StaticMatrix, SArray, SVector, SMatrix, SOneTo

using PrettyPrinting
const Pretty = PrettyPrinting

import FillArrays
using Static
using Static: StaticInteger
using FunctionChains
using PropertyFunctions: PropSelFunction

using StaticThings:
    SizeLike, StaticSizeLike, OneToLike, StaticOneToLike, IntegerLike,
    asaxes, asnonstatic, canonical_size, size_tuple,
    maybestatic_eachindex, maybestatic_length, maybestatic_size,
    maybestatic_first, maybestatic_last, maybestatic_view,
    maybestatic_oneto, maybestatic_fill, maybestatic_reshape,
    axes2size, size2length, maybestatic_split, staticarray_type,
    tuple_type_all, tuple_type_any, tuple_type_mapreduce,
    sum_leading_dims, drop_leading_dims, merge_leading_dims, all_leading_dims

import HeterogeneousComputing
using HeterogeneousComputing: real_numtype
using HeterogeneousComputing:
    GenContext, get_rng, get_precision, get_compute_unit, is_host_unit, allocate_array, fill_array

using ArraysOfArrays:
    ArrayOfSimilarArrays, VectorOfArrays, VectorOfVectors, VectorOfSimilarArrays,
    VectorOfSimilarVectors, flatview, fused, stacked, sliced, getsplitmode,
    is_memordered_splitmode, AbstractSplitMode, AbstractPartMode, UnknownSplitMode, NonSplitMode

using OneTwoMany: firstarg, secondarg

using StructArrays: StructArray
import StructArrays
import Adapt

export gentype

export AbstractMeasure

import IfElse: ifelse
export logdensity_def
export basemeasure
export basekernel
export productmeasure

export insupportof
export getdof
export transport_to

include("insupportof.jl")

abstract type AbstractMeasure end

AbstractMeasure(m::AbstractMeasure) = m

"""
    asmeasure(m)

Turns a measure-like object `m` into an `AbstractMeasure`.

Calls `convert(AbstractMeasure, m)` by default
"""
function asmeasure end

@inline asmeasure(m::AbstractMeasure) = m
asmeasure(m) = convert(AbstractMeasure, m)
export asmeasure

"""
    struct AsMeasure{T}

Wraps a measure-like object into an `AbstractMeasure`.

Constructor:

```
AsMeasure{T}(obj::T)
```

User code should not create instances of `AsMeasure` directly, but should
call `asmeasure(obj)` instead.
"""
struct AsMeasure{T} <: AbstractMeasure
    obj::T

    AsMeasure{T}(obj::T) where {T} = new(obj)
end

# Struct arrays of wrapped objects rebuild elements via ConstructionBase:
ConstructionBase.constructorof(::Type{<:AsMeasure}) = _asmeasure
_asmeasure(obj) = AsMeasure{typeof(obj)}(obj)

Base.:(==)(a::AsMeasure, b::AsMeasure) = a.obj == b.obj
Base.hash(a::AsMeasure, h::UInt) = hash(a.obj, hash(:AsMeasure, h))
Base.isapprox(a::AsMeasure, b::AsMeasure; kwargs...) = isapprox(a.obj, b.obj; kwargs...)

function Pretty.quoteof(d::M) where {M<:AbstractMeasure}
    the_names = fieldnames(typeof(d))
    :($M($([getfield(d, n) for n in the_names]...)))
end

@inline DensityKind(::AbstractMeasure) = HasDensity()

Broadcast.broadcastable(m::AbstractMeasure) = Ref(m)

gentype(μ::AbstractMeasure) = typeof(testvalue(μ))

# gentype(μ::AbstractMeasure) = gentype(basemeasure(μ))

using LogExpFunctions: logsumexp, logistic, logit

@deprecate instance_type(x) Core.Typeof(x) false


"""
`logdensity_def` is the standard way to define a log-density for a new measure.
Note that this definition does not include checking for membership in the
support; this is instead checked using `insupportof`. `logdensity_def` is
a low-level function, and should typically not be called directly. See
`logdensityof` for more information and other alternatives.

---

    logdensity_def(m, x)

Compute the log-density of the measure m at the point `x`, relative to
`basemeasure(m)`, and assuming `insupportof(m, x)`.

---

    logdensity_def(m1, m2, x)

Compute the log-density of `m1` relative to `m2` at the point `x`, assuming
`insupportof(m1, x)` and `insupportof(m2, x)`.

The generic implementation descends the base measure chains of both
measures in lockstep, so it terminates at the first pair of base measures
for which a specialized relative density is available (in particular at
pairs of identical primitive measures) and any shared chain suffix cancels
symbolically. To provide specialized relative densities for pairs of
measure types, add methods to [`MeasureBase.logdensity_rel_def`](@ref),
not to `logdensity_def` itself.
"""
function logdensity_def end

using Compat

using IrrationalConstants
using IrrationalConstants: loghalf

include("collection_utils.jl")
include("smf.jl")
include("mspace.jl")
include("getdof.jl")
include("standard/stdmeasure.jl")
include("transport.jl")
include("rand.jl")
include("proxies.jl")
include("parameterized.jl")
include("domains.jl")
include("primitive.jl")
include("utils.jl")
include("special-derivatives.jl")
include("mass-interface.jl")

include("density.jl")
include("density-core.jl")
include("density-batched.jl")

include("primitives/counting.jl")
include("primitives/lebesgue.jl")
include("primitives/dirac.jl")
include("primitives/trivial.jl")

include("combinators/transformedmeasure.jl")
include("combinators/reshape.jl")
include("combinators/weighted.jl")
include("combinators/normalized.jl")
include("combinators/superpose.jl")
include("combinators/product.jl")
include("combinators/power.jl")
include("transport-batched.jl")
include("combinators/spikemixture.jl")
include("combinators/likelihood.jl")
include("combinators/restricted.jl")
include("combinators/smart-constructors.jl")
include("combinators/conditional.jl")
include("combinators/implicitlymapped.jl")

include("standard/stduniform.jl")
include("standard/stdexponential.jl")
include("standard/stdlogistic.jl")
include("standard/stdnormal.jl")
include("standard/stdconvert.jl")
include("standard/stdtraits.jl")
include("combinators/combined.jl")
include("combinators/bind.jl")
include("combinators/half.jl")
include("combinators/truncated.jl")

#include("implicitmaps.jl")

include("measure_operators.jl")

using .MeasureOperators
export ⋄, ⊙, ▷, ⊗, ∫, ∫exp, 𝒹, log𝒹

include("interface.jl")

using .Interface

# Public but not exported: extension points and internals for measure
# implementations.
@compat public MeasureOperators,
    # Densities:
    logdensityof_impl, batched_logdensityof_impl, batched_logdensity_def,
    logdensityof_with_rest, batched_logdensityof_with_rest,
    logdensity_rel_def, as_integrand, as_integrand_exp,
    # Random variates:
    GenContext, rand_impl, batched_rand_impl,
    # Transport:
    transport_def, batched_transport_def,
    transport_to_std, transport_from_std,
    transport_to_std_with_rest, transport_from_std_with_rest,
    batched_transport_to_std, batched_transport_from_std,
    batched_transport_to_std_with_rest, batched_transport_from_std_with_rest,
    preferred_stdmeasure, promote_stdmeasure, StdMeasure, AnyStdMeasure,
    NoStdTransport, stdconvert, TransportFunction,
    # Measurable spaces, degrees of freedom and argument checks:
    mspace_ndims, mspace_flatsize, mspace_elsize, NoMSpaceElementSize, fixed_stream_size,
    checked_arg, NoDOF, NoFastDOF, NoFastInsupport,
    # Masses:
    AbstractUnknownMass, UnknownMass, UnknownFiniteMass,
    # Measure types and accessors:
    pwr_base, pwr_axes, pwr_size, DensityMeasure, AsMeasure, StdSimplex,
    transportmeasure, tpmeasure_split_combined,
    gettransform, TransformVolCorr, NoVolCorr, WithVolCorr

end # module MeasureBase
