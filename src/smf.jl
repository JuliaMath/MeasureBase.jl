@doc raw"""
    smf(μ, x::Real) ::Real

Compute the _Stieltjes measure function (SMF)_ of the measure `μ` at the point
`x`.

The SMF is the measure-theoretic generalization of the _cumulative distribution
function (CDF)_ from probability theory. An SMF `F(x) = smf(μ, x)` must have the
following properties:

1. F is _nondecreasing_
2. F is _right-continuous_: `F(x)` should be the same as `lim_{δ→0} F(x + |δ|)`.
3. μ((a,b]) = F(b) - F(a)

Note that unlike the CDF, an SMF is only determined up to addition by a
constant. For many applications, this leads to a need to evaluate an SMF at -∞.
It's therefore important that `smf(μ, -Inf)` be fast. In practice, this will
usually be called as `smf(μ, static(-Inf))`. It's then easy to ensure speed and
avoid complex control flow by adding a method `smf(μ::M, ::StaticFloat64{-Inf})`.

Users who pronounce `sinh` as "sinch" are advised to pronounce `smf` as "smurf".
"""
function smf end

export smf

function smfinv end

export smfinv

struct NoSMF end

struct NoSMFInverse end

smf(μ, x) = NoSMF()

smfinv(μ, p) = NoSMFInverse()