"""
    isprimtype(μ)

Most measures are defined in terms of other measures, for example using a
density or a pushforward. Those that are not are considered (in this library,
it's not a general measure theory thing) to be _primitive_. The canonical
example of a primitive measure is `Lebesgue(X)` for some `X`.

The default method is
    isprimtype(μ) = false

So when adding a new primitive measure, it's necessary to add a method for its type
that returns `true`.
"""
function isprimtype end

@traitdef IsPrimType{X}
@traitdef IsPrimType{X} <- isprimtype(X)
isprimtype(X) = false # default

export basemeasure

"""
    basemeasure(μ)

Many measures are defined in terms of a logdensity relative to some base
measure. This makes it important to be able to find that base measure.

For measures not defined in this way, we'll typically have `basemeasure(μ) == μ`.
"""
function basemeasure end

basemeasure(μ::M) where {IsPrimType{M}} = μ

include("primitives/trivial.jl")
include("primitives/lebesgue.jl")
include("primitives/counting.jl")
include("primitives/dirac.jl")

@traitdef IsRepType{X}

@traitimpl IsRepType{X} <- isreptype(X)
isreptype(X) = false # set default

# Every primitive measure is also a representative
isreptype(μ::M) where {IsPrimType{M}} = true
