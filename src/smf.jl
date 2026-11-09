@doc raw"""
    smf(μ, x::Real) ::Real
    smf(μ, a::Real, x::Real) ::Real

Compute the _Stieltjes measure function (SMF)_ of the measure `μ` at the point
`x`, resp. the mass `μ((a, x])` of an interval.

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

`smf(μ, a, x)` is the signed difference `smf(μ, x) - smf(μ, a)`, the mass
`μ((a, x])` for `a <= x` and `-μ((x, a])` for `x < a`, but computed
without the cancellation of the difference where `x` is close to `a`. This
is the extension point for accurate masses of short intervals, which
truncations and [`Half`](@ref) build on. The default is the difference of
`smf` (of `csmf` if `a` is in the upper tail, see [`csmf`](@ref)).
Standard measures specialize it.

Users who pronounce `sinh` as "sinch" are advised to pronounce `smf` as "smurf".
"""
function smf end

export smf

"""
    invsmf(μ, p)
    invsmf(μ, a, p)

The inverse of the Stieltjes measure function [`smf`](@ref) of the measure
`μ`, `smf(μ, invsmf(μ, p)) ≈ p` for measures without atoms.

`invsmf(μ, a, p)` is the inverse of the interval form `smf(μ, a, x)` in
`x`, accurate for `x` close to `a`: the `x` above `a` with `μ((a, x]) = p`
for `p >= 0`, the `x` below `a` with `μ((x, a]) = -p` for `p < 0`.
`NaN` for `p` outside of `[-μ((-∞, a]), μ((a, ∞))]`.
"""
function invsmf end

export invsmf

@doc raw"""
    csmf(μ, x::Real)::Real

The complementary Stieltjes measure function of the measure `μ` of finite
mass, the measure `μ((x, ∞))` of the interval above `x`, the
generalization of the complementary cdf (ccdf).

`csmf(μ, x)` is `massof(μ) - smf(μ, x)` for Stieltjes measure functions
that vanish at `-∞`, but computed accurately in the upper tail, where
`smf(μ, x)` rounds to the mass of `μ`. `μ((a, b]) = csmf(μ, a) - csmf(μ, b)`
is accurate for intervals in the upper tail, `smf(μ, b) - smf(μ, a)` for
others.

See also [`invcsmf`](@ref), [`logsmf`](@ref) and [`logcsmf`](@ref).
"""
function csmf end

export csmf

"""
    invcsmf(μ, p)

The inverse of the complementary Stieltjes measure function [`csmf`](@ref)
of the measure `μ`, `csmf(μ, invcsmf(μ, p)) ≈ p` for measures without atoms,
accurate for small `p` (the upper tail).
"""
function invcsmf end

export invcsmf

"""
    logsmf(μ, x)

The logarithm of the Stieltjes measure function [`smf`](@ref), accurate
also where `smf(μ, x)` underflows (the far lower tail). Defaults to
`log(smf(μ, x))`, see also [`invlogsmf`](@ref).
"""
function logsmf end

export logsmf

"""
    logcsmf(μ, x)

The logarithm of the complementary Stieltjes measure function
[`csmf`](@ref), accurate also where `csmf(μ, x)` underflows (the far upper
tail). Defaults to `log(csmf(μ, x))`, see also [`invlogcsmf`](@ref).
"""
function logcsmf end

export logcsmf

"""
    invlogsmf(μ, logp)

The inverse of [`logsmf`](@ref), defaults to `invsmf(μ, exp(logp))`.
"""
function invlogsmf end

export invlogsmf

"""
    invlogcsmf(μ, logp)

The inverse of [`logcsmf`](@ref), defaults to `invcsmf(μ, exp(logp))`.
"""
function invlogcsmf end

export invlogcsmf

struct NoSMF end

struct NoSMFInverse end

smf(μ, x) = NoSMF()
invsmf(μ, p) = NoSMFInverse()
csmf(μ, x) = NoSMF()
invcsmf(μ, p) = NoSMFInverse()

# Intervals starting in the upper tail (where `csmf(μ, a) < smf(μ, a)`) have
# their mass resolved by the complementary Stieltjes measure function:
@inline _upper_side(F, G) = G < F
@inline _upper_side(::Any, ::NoSMF) = static(false)

function smf(μ, a, x)
    F_a, G_a = smf(μ, a), csmf(μ, a)
    _select(_upper_side(F_a, G_a), () -> G_a - csmf(μ, x), () -> _smf_diff(smf(μ, x), F_a))
end

@inline _smf_diff(F_x, F_a) = F_x - F_a
@inline _smf_diff(::NoSMF, ::NoSMF) = NoSMF()

function invsmf(μ, a, p)
    F_a, G_a = smf(μ, a), csmf(μ, a)
    _select(_upper_side(F_a, G_a), () -> invcsmf(μ, G_a - p), () -> _invsmf_above(μ, F_a, p))
end

@inline _invsmf_above(μ, F_a, p) = invsmf(μ, F_a + p)
@inline _invsmf_above(μ, ::NoSMF, p) = NoSMFInverse()

logsmf(μ, x) = _log_smf(smf(μ, x))
logcsmf(μ, x) = _log_smf(csmf(μ, x))
@inline _log_smf(F) = log(F)
@inline _log_smf(F::NoSMF) = F

invlogsmf(μ, logp) = invsmf(μ, exp(logp))
invlogcsmf(μ, logp) = invcsmf(μ, exp(logp))

smf(μ) = Base.Fix1(smf, μ)
invsmf(μ) = Base.Fix1(invsmf, μ)
