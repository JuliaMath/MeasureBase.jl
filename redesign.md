# MeasureBase redesign notes (branch `major-upgrade`)

Working notes on the approach behind this branch, for reviews and for
guiding the next steps. Kept up to date while the branch evolves, to be
removed before the merge.

## Goals

- A breaking release that runs on GPUs (CUDA, JLArrays) and under
  Reactant, with batching built into the foundation.
- One implementation per measure type for densities, transports and
  random variates, so scalar and batched paths can't drift apart.
- Composable, structural solutions: powers, products, combinations,
  binds and pushforwards implement their behavior once in terms of their
  components. No shape inference, no per-call DOF sums, no function
  traits in the core.

## Design philosophy

- **Batched first.** A single variate is a batch with zero batch
  dimensions. Every kernel handles both, the point API is the batched API
  at zero batch dimensions. Static arrays keep the scalar path
  allocation-free.
- **Ranks, not sizes.** A flat batch is an array
  `(variate dims..., batch dims...)`. Kernels only need the variate rank
  of their measure; sizes are optional declarations for validation,
  stream consumption and static fast paths, never for routing. Unknown
  sizes are safe.
- **Standard measures as pivots.** Transport goes through a standard
  measure type chosen by promoting the measures' preferences. Measure
  types implement transport to and from their preferred standard measure
  only.
- **Branch-free, device-friendly kernels.** Broadcasts, reductions and
  masks instead of branches; host loops only where documented.
- **Entry points normalize layouts.** Users pass nested arrays, tuples of
  batches, struct arrays or flat arrays; the kernels only see flat
  batches.

## Concepts

**Variate rank.** `mspace_ndims(::Type{M})`, 0 for scalar variates.
Derived from `mspace_flatsize(::Type{M})` where known, declared by
array-variate leaves, derived by structural measures. Without a rank the
batched defaults throw with a message naming the declaration; point
kernels keep working.

**Flat batches.** Kernels return arrays over the batch dimensions, a
number for a single variate, possibly lazily. ArraysOfArrays containers
are fused into their flat storage at the entry points, ragged
containers are evaluated variate by variate. Tuple and named tuple
variates batch as tuples of batches; struct arrays and arrays of tuples
are accepted, their flat storage is the tuple of component storages.
`rand(Pt^n)` of a tuple product is a struct array.

**Streams.** Variates of `mcombine(vcat, ...)`, binds and tuple products
inside such streams are flat vectors consumed with the with-rest
protocol. Point forms return `(result, x_μ, x_rest)` (binds need the
consumed variate), batched forms take streams `(rows, batch dims...)`
and a multiplicity `sz::SizeLike` of variates per stream and return
`(result, rest)`. Powers pass their size as multiplicity to their base;
combined measures and tuple products split rows by their fixed stream
lengths. `fixed_stream_size(::Type{M})` decides whether a batch of
streams is consumed in fused operations or stream by stream by the
outermost combinator (binds never fuse). Scalar leaves consume
`(1, sz..., batch dims...)` and drop the leading dimension. Nested
element variates in vcat streams are flattened.

**Transport.** `transport_to(ν, μ)` gives a `TransportFunction`; `f(x)`
transports a variate, `f.(X)` a whole batch, `transports_with_ladj(f, X)`
also returns the log-abs-det Jacobians over the batch (density ratios
via `logdensities`). Results never alias the input. Leaves implement
`transport_to_std`/`transport_from_std` for their preferred standard
measure (`preferred_stdmeasure`, `promote_stdmeasure`, `AnyStdMeasure`,
`NoStdTransport`), array-variate leaves also the `batched_` forms.
`stdconvert` converts between standard measures in log form,
`transport_def` may be specialized for direct pairs. Transports map
normalized laws: weighted measures transport like their base, and the
log-abs-det Jacobian is the density ratio of the normalized measures, with
weights stripped structurally (through weighted measures, tuple products,
powers, combinations, pushforwards and binds) and the remaining leaves
normalized by their mass. Leaves of unknown mass make
`with_logabsdet_jacobian` throw. Standard streams are
`(dof, batch dims...)`; the from-side default consumes `fast_dof(μ)`
entries per variate, composed measures implement the with-rest forms.
`fast_dof` is static or constant-time (standard measures, Lebesgue, Dirac,
wrapped univariate and Dirichlet distributions, powers, tuple products,
arrays of single-DOF marginals, combinations), `NoFastDOF` otherwise
(default: the fast DOF of the base measure, static `getdof` of
primitives). Transports and draws split streams with it, without it
composed measures take the with-rest path and leaves fall back to
`getdof` once per call. Neither is used inside kernels.

**Static variates.** Sizes are `StaticThings.SizeLike` throughout, the
with-rest multiplicity and the batch size of `rand` included, so a
statically sized power stays static through streams. On the CPU, bulk
draws of a fully static size are static arrays and the stream splits and
reshapes keep them static: `rand`, `transport_to` and `logdensityof` of a
measure whose variate sizes are all static are type stable and allocation
free, for powers, tuple and named tuple products and binds over them.
Sizes that depend on variate values stay dynamic, and other compute units
allocate their own arrays.

**Random variates.** `rand(ctx::GenContext, μ)` with RNG, precision and
compute unit (`rand(μ)`, `rand(rng, μ)`, `rand(T, μ)` are wrappers).
`batched_rand_impl(ctx, μ, sz::SizeLike)` returns a flat batch, a single
variate for `sz == ()`; `rand_impl` defaults to it. Defaults draw
standard variates in bulk on the compute unit and transport them, or
generate variate by variate without a standard transport. Draws go
through HeterogeneousComputing (0.3): `rand`/`randn`/`randexp(ctx, dims)`
for bulk draws and `rand(ctx)` etc. for single ones, on any compute unit
and RNG (exponential draws from uniform ones where RNGs lack them, scalar
draws from device RNGs), `fill_array` for constants. MeasureBase only
adds static arrays for static sizes on host units and constant fills for
`ConstantRNG` test values, which HC's array draws don't support. Under
Reactant the context, with a `ReactantRNG`, is an argument of the
compiled function.

**Array products.** `productmeasure(::AbstractArray)` stores isbits
parameterized marginals as `StructArrays` (nested parameter structs
unwrapped; numbers, arrays, tuples, strings, symbols and function objects
stay opaque columns), in one place: `_marginal_storage`. Fused kernels
broadcast over the leaf columns and rebuild marginals via
`ConstructionBase.constructorof`, which works on CUDA and under Reactant.
Fusion needs concrete scalar-variate marginals (one DOF for transport);
other array products loop over host-resident marginals. Measures holding
arrays have `Adapt` rules.

**Parameter batches.** Arrays of measures are evaluated with
`logdensities(ms, X)` (batch dimensions of `X` broadcast against `ms`)
and by broadcasting `logdensityof`/`densityof` over them, which share the
implementation: scalar-variate measures in one broadcast over the leaf
columns of struct arrays (`_marginal_broadcast`), powers of one shape via
their bases over the power dimensions, other array variates measure by
measure. `logdensityof.(Ref(ℒ), Θ)` builds the kernel's measures by
broadcasting over `Θ` (struct arrays stay struct arrays), scalar
observations broadcast over the parameter columns directly (traced too).
`DensityMeasure` batches are base batch plus integrand batch, the
integrand evaluated on the in-support subset on the host (everywhere,
masked, on devices and when traced).

**Pushforwards.** `pushfwd(f, μ)` never evaluates `f` at construction.
The output size is declared (`pushfwd(f, μ; varsize)`), the origin's for
`Base.BroadcastFunction`s, or inferred from the type `f` returns for
numbers or arrays of the origin's rank (static sizes for static origins,
`Core.Compiler.return_type`, used only for these declarations); otherwise
it is learned from `f(testvalue(μ))` on demand (not cached, point
densities don't need it). Batched
application needs no traits: elementwise for `Base.BroadcastFunction`
(fused density kernels), column batches for AffineMaps types (weak
dependency, `MeasureBaseAffineMapsExt`), a host loop otherwise. Measures
with counting root measures push their root measure forward (counting
measure on the image) instead of adapting it with a volume element, also
under `AdaptRootMeasure()`, and so do measures on the simplex (its
coordinate measure is a chart, maps of the ambient space have no volume
element for it); mixed roots (superpositions of continuous and atomic
measures) still get the volume element everywhere. Pushforwards by
transports adapt the root measure to the one of the transport's target
(the volume element of a transport is the density ratio relative to the
roots of its source and target) and don't compose with other
pushforwards.

## Extension points

| Aspect | Scalar-variate leaf | Array-variate leaf | Composed measure |
|---|---|---|---|
| Density | `logdensity_def` | `mspace_ndims`, `batched_logdensityof_impl` (+`_def`) | kernels via components, with-rest forms |
| Transport | `preferred_stdmeasure`, `transport_to_std`, `transport_from_std` | + `batched_transport_to_std`, `batched_transport_from_std` | with-rest forms, `fixed_stream_size` |
| Random | `batched_rand_impl` (default via transport), optionally `rand_impl` | `batched_rand_impl` | derived |
| Declarations | none | `mspace_ndims`, optionally `mspace_flatsize`, `getdof` | derived |

## Layouts per combinator

- Powers: `(base dims..., power dims..., batch dims...)`, innermost base
  first; standard variates are the flat vector of the base's. Results
  are nested views over the flat storage. Powers of tuple products treat
  numeric flat variates as streams and accept tuples of batches.
- Array products: `(marginal dims..., product dims..., batch dims...)`.
- Tuple products: tuples of batches; marginal by marginal in streams.
- Combined `vcat`: streams; `merge`: merged named tuples. Binds: single
  streams, value-dependent sizes.
- Weighted, restricted, density measures, Half: forward plus weights or
  masks. Superpositions and spike mixtures: one batch per component,
  masks aligned with the variate dimensions. Dirac: constant batches.
- Distributions extension: the main univariate families (Normal,
  Uniform, Exponential, Logistic, Cauchy, Laplace, LogNormal, Weibull,
  Gamma, Beta, Poisson, Bernoulli) get branch-free density formulas and
  transports to the standard measure matching their shape (`families.jl`),
  other univariate distributions go via `StdLogistic` (log-cdf and
  quantile). `MvNormal` works on column batches through its Cholesky
  factor, `Dirichlet` through stick-breaking Beta transports in FlatPPL's
  orientation (break fractions `Beta(α_i, Σ_{j>i} α_j)`, the complement
  of BAT 5.x's coordinates) with `cumsum` and `cumprod` along the variate
  dimension (`multivariate.jl`), with
  densities relative to `Lebesgue(StdSimplex(n))`, the coordinate measure
  of the simplex (as in FlatPPL), its own base and root measure. Incomplete
  gamma and beta functions go through the `_gamma_cdf`/`_beta_cdf` hooks
  (and quantiles), whose ForwardDiff and ChainRules derivatives come from
  the densities and, for the shape parameters, from
  `_gamma_cdf_dshape`/`_beta_cdf_dshape` (term-wise differentiated series
  and continued fractions, Moore 1982 and Boik & Robison-Cox 1998; host
  only, like SpecialFunctions). Draws use the Distributions samplers on
  host compute units (`is_host_unit`) and the standard transports on
  others, and on the host too where the parameter precision (`partype`) is
  below the context precision. Distributions without standard transports
  convert their draws afterwards (the context then sets the type, not the
  resolution). Remaining array-variate distributions batch via
  `logpdf(d, X)` on the host.

## Changes relative to `main`

After the initial rework of domains, densities, likelihoods and binds
(with the Distributions and autodiff extensions), the branch went on with
the variate size contract and `preferred_stdmeasure`, densities over flat
storage, branch-free evaluation, Reactant smoke tests, structural batched
kernels, transport rebuilt on standard measures, batched transport and the
broadcast hook, rand via `GenContext`, then the batched-first redesign
(density core, struct array products, transport, random variates,
structured batches), the review fixes, and device-friendly kernels for
wrapped Distributions.

Removed: `transport_origin`/`to_origin`/`from_origin` and the origin
machinery, `NoTransport`, `transport_to_mvstd`, per-measure
`Base.rand(rng, T, μ)` methods, the `_trafo_cdf`/`_trafo_quantile` hooks,
and the exports `kleisli`, `AbstractTransitionKernel`, `↣`, `rebase`,
`commonbase`, `likelihood_ratio`, `log_likelihood_ratio`, `↑`,
`ZeroSet`, `CodimOne`, `Simplex`, `ℝ₊` and `𝕀` (with
`PointwiseProductMeasure` and power-weighted measures).

API: the operators live in `MeasureOperators` and are re-exported; `⊙`
is now the pullback (it was the pointwise product). Extension points and
the internals BAT uses are declared `public` (`@compat public`, effective
on Julia 1.11+). `mspace_elsize` (size of
array variates, `NoMSpaceElementSize` for structured ones) is only used
for with-rest stream consumption and `mreshape`, and is to be unexported
(and declared public). Kernels use `mspace_ndims` and `mspace_flatsize`
(size of the flat storage, which differs for nested powers). `some_dof`
(DOF at an unspecified point) sizes standard transport partners of
measures without a fixed DOF. The release notes go into the pull request
description, the concepts and extension points are in docs/src/concepts.md.

Behavior changes for the release notes: univariate Distributions pivot on
`StdLogistic`; nested powers and array products return ArraysOfArrays
views, powers of tuple products struct arrays; `Half` transports via
`StdUniform` (tails limited); `mcombine(vcat, ::Product, ::Product)`
merges only concrete homogeneous marginals; `rand(rng, Int, μ)`
unsupported; isbits marginal arrays become struct arrays (`rand` of such
products gives plain arrays); vcat-combined and bind variates are flat;
`transport_to(StdUniform, m)` with binds inside errors (use
`StdUniform()^n`); batched kernels of rank-less measures throw.

## Verification

Full suite (Aqua, extensions, doctests) on CPU with JLArrays cases and
a FixedSizeArrays transparency check (fixed-size inputs give fixed-size
outputs; FixedSizeArrays stays a test dependency, allocating fixed-size
variates by default is a HeterogeneousComputing decision for later).
`test/static_variates.jl` checks that statically sized measures stay type
stable and allocation free.
`test/test_reactant.jl` runs as part of the suite on 64-bit Linux and
macOS with stable Julia, adding Reactant on demand as MGVI does (backend
via `MEASUREBASE_REACTANT_BACKEND`). `test/cuda` is opt-in. Both run
locally on the GB10, green at HEAD except one expected-broken CUDA case
(AffineMaps Jacobian on device).

## Known gaps and open decisions

- Upstream: AffineMaps lacks `Adapt` rules and device/traced Jacobians;
  Distributions isn't device-aware and its parameter structs (PDMats,
  `Dirichlet`) can't hold traced arrays, so under Reactant distribution
  parameters stay constants; SpecialFunctions' incomplete beta and gamma
  functions don't compile for GPUs or Reactant, so Beta and Dirichlet
  transports and draws (Gamma under Reactant too) stay on the CPU;
  ChangesOfVariables has no rules for `Base.Fix1`/`Fix2` arithmetic;
  HeterogeneousComputing has no splittable per-stream RNGs yet and
  doesn't allocate FixedSizeArrays; JLArrays has no triangular solves;
  Reactant rejects traced `VectorOfArrays` and empty batches.
- Out-of-support convention (decided 2026-09-19): structural errors (wrong
  rank, size, container or element kind) throw an `ArgumentError` from
  `checked_arg` at the entry points; variates of the right shape in the
  variate space of the measure never throw: densities are `-Inf` outside
  the support (including non-integers for counting-based measures and
  infinite values), transports are `NaN` outside the support of the
  source, `NaN` inputs give `NaN` or `-Inf`, relative densities keep
  `+Inf`/`-Inf`/`NaN`. Pushforwards are only defined on the image of their
  map, the inverse may throw outside of it (see below). Finite inputs give
  finite transports: uniform inputs are clamped into the open unit interval
  before quantiles (the endpoints stand for their nearest interior grid
  points, a null-set convention like BAT's `[eps, 1 - eps]` clamping in
  practice), tail probabilities of log-space conversions floor at the
  smallest positive float.
  Kernels must not throw outside the support, since the masks evaluate
  both branches (`abs`, `clamp`, `min` guards instead of `NaNMath`, which
  isn't device-compatible). Downstream checks such as BAT's
  `checked_logdensityof` stay downstream.
- BAT's boundary tweaks: adopted are infinite variates outside the
  support of continuous wrapped distributions, the clamping of uniform
  inputs into the open unit interval (for inputs inside `[0, 1]`, outside
  stays `NaN`), `-Inf + Inf` giving `-Inf` in pushforward densities and a zero
  Jacobian term where both densities vanish. Not adopted: replacing finite
  densities with infinite Jacobian terms by `-1e38`, and re-evaluating
  densities an `eps` inside the support where Distributions returns `NaN`
  (the family kernels are exact there).
  Quantile results within `4 eps` of the support edges snap to the edges
  only on the generic logistic path of wrapped distributions.
- Base measures of discrete measures (review fix, 2026-09-28): a `Dirac`
  is its own base measure and `Counting(S)` keeps `S` as its base, so the
  base-measure mask of superpositions sees the atoms and a superposition
  of a continuous and a discrete component has densities relative to
  Lebesgue measure plus counting measure on the atoms. `CountingBase()` on
  all of the reals is not σ-finite and is never used as a base of a
  superposition component.
- Bind sampling and masses (decided 2026-10-22): `rand` of a bind samples
  ancestrally from the normalized laws of its components (the primary
  variate from the normalized primary measure, the secondary one from the
  normalized kernel measure), the normalized law of the bind only for
  kernels of constant mass. The mass of a bind integrates the kernel
  masses over the primary measure, known for weighted and superposed point
  masses, unknown otherwise. Projections (`mcombine` and `mbind` with
  `firstarg`/`secondarg`) weight the kept measure with the mass of the
  dropped one (`massof(α) * massof(β)`), pointwise for binds; where that
  mass is unknown, so are the mass and the density of the projection.
- Pushforwards through non-surjective maps (decided 2026-10-23):
  pushforwards are only defined on the image of their map, the inverse may
  throw outside of it (`logdensityof(pushfwd(exp, StdNormal()), -1.0)`
  throws a `DomainError` from `log`). Documented, not masked.
- Empty batches of array products over heterogeneous marginals fail in
  the reduction (`zero(Any)`), needs an `init` keyword for
  `sum_leading_dims` in StaticThings.
- `convert_realtype` turns Poisson/Bernoulli draws into floats while
  `Counting(Int)` checks types, so `logdensity_rel(Poisson, Counting(Int),
  2.0)` is `+Inf` but finite relative to `Counting(ℤ)`. Either keep
  integer draws or make `Counting(Int)` check values.
- Decisions pending: `Half` general design (law of `abs(x)`, density
  `logaddexp(ℓ(x), ℓ(-x))`, quantile by root finding unless symmetric,
  symmetry declared by the Distributions extension) versus the current
  symmetric-only `Half`, `Half` tails via log-ccdf, the
  `smart-constructors.jl` review (location-scale arrays as affine
  pushforwards of powers), `_static_ndims` type-first vs. instance-first.
- `normalize` (decided 2026-10-23, review D7): `normalize(μ)` is a
  `NormalizedMeasure(base, logmass)`, a weighted measure with the weight
  `-logmass`, or `μ` itself if `massof(μ)` is statically one. It strips
  weights and is idempotent, zero and infinite masses throw. Of unknown
  mass (`logmass` is an `AbstractUnknownMass`) it has no density
  (`ArgumentError`), but samples, transports and has the support of its
  base. `massof` is type stable and known for `Lebesgue(a..b)`, `Counting`
  of bounded integers, ranges and sets, weighted measures, products,
  powers and superpositions, restrictions stay unknown. Transport Jacobians
  normalize the leaves of the structural weight stripping, see Transport.
- Complementary Stieltjes measure functions (decided 2026-10-23): `csmf`
  (`μ((x, ∞))`, ccdf-like, accurate in the upper tail), `invcsmf`, and the
  log forms `logsmf`, `logcsmf`, `invlogsmf`, `invlogcsmf` (defaults via
  `log`/`exp`, stable forms for the standard measures and wrapped
  distributions). `truncate` and interval masses use `csmf` where
  `csmf(μ, a) < smf(μ, b)` (branch-free for traced bounds).
- Constant batches (2026-10-23): batches of point masses with numeric atoms
  are `Fill`s on the host, from `rand` and transports alike, dense arrays
  on devices.
- Names clashing with Distributions (decided 2026-10-23): `insupport` and
  `params` are renamed to `insupportof` and `paramsof`, the old names stay
  exported as deprecations for one release.
- Polish before merge: final docs and release notes pass, history curation,
  version bump, remove this file.
