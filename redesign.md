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
and a multiplicity `sz::Dims` of variates per stream and return
`(result, rest)`. Powers pass their size as multiplicity to their base;
combined measures and tuple products split rows by their fixed stream
lengths. `fixed_stream_size(::Type{M})` decides whether a batch of
streams is consumed in fused operations or stream by stream by the
outermost combinator (binds never fuse). Scalar leaves consume
`(1, sz..., batch dims...)` and drop the leading dimension. Nested
element variates in vcat streams are flattened.

**Transport.** `transport_to(ν, μ)` gives a `TransportFunction`; `f(x)`
transports a variate, `f.(X)` a whole batch. Leaves implement
`transport_to_std`/`transport_from_std` for their preferred standard
measure (`preferred_stdmeasure`, `promote_stdmeasure`, `AnyStdMeasure`,
`NoStdTransport`), array-variate leaves also the `batched_` forms.
`stdconvert` converts between standard measures in log form,
`transport_def` may be specialized for direct pairs. Standard streams are
`(dof, batch dims...)`; the from-side default consumes `fast_dof(μ)`
entries per variate, composed measures implement the with-rest forms.
`getdof`/`fast_dof` are declaration-derived, used at construction time
and for chunking, never inside kernels.

**Random variates.** `rand(ctx::GenContext, μ)` with RNG, precision and
compute unit (`rand(μ)`, `rand(rng, μ)`, `rand(T, μ)` are wrappers).
`batched_rand_impl(ctx, μ, sz::Dims)` returns a flat batch, a single
variate for `sz == ()`; `rand_impl` defaults to it. Defaults draw
standard variates in bulk on the compute unit and transport them, or
generate variate by variate without a standard transport.

**Array products.** `productmeasure(::AbstractArray)` stores isbits
parameterized marginals as `StructArrays` (nested parameter structs
unwrapped; numbers, arrays, tuples, strings, symbols and function objects
stay opaque columns), in one place: `_marginal_storage`. Fused kernels
broadcast over the leaf columns and rebuild marginals via
`ConstructionBase.constructorof`, which works on CUDA and under Reactant.
Fusion needs concrete scalar-variate marginals (one DOF for transport);
other array products loop over host-resident marginals. Measures holding
arrays have `Adapt` rules.

**Pushforwards.** `pushfwd(f, μ)` learns its output size once at
construction from a test value when the origin has a size. Batched
application needs no traits: elementwise for `Base.BroadcastFunction`
(fused density kernels), column batches for AffineMaps types (weak
dependency, `MeasureBaseAffineMapsExt`), a host loop otherwise.

## Extension points

| Aspect | Scalar-variate leaf | Array-variate leaf | Composed measure |
|---|---|---|---|
| Density | `logdensity_def` | `mspace_ndims`, `batched_logdensityof_impl` (+`_def`) | kernels via components, with-rest forms |
| Transport | `preferred_stdmeasure`, `transport_to_std`, `transport_from_std` | + `batched_transport_to_std`, `batched_transport_from_std` | with-rest forms, `fixed_stream_size` |
| Random | `batched_rand_impl` (default via transport) | `batched_rand_impl` | derived |
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
- Distributions extension: univariate via `StdLogistic` (log-cdf and
  quantile), location-scale families via their affine map, `MvNormal` via
  Cholesky factors, array-variate batches via `logpdf(d, X)`.

## Changes relative to `master`

Twenty commits since `c773afe`: variate size contract and
`preferred_stdmeasure`, densities over flat storage, branch-free
evaluation, Reactant smoke tests, structural batched kernels, transport
rebuilt on standard measures, batched transport and the broadcast hook,
rand via `GenContext`, then the batched-first redesign (density core,
struct array products, transport, random variates, structured batches)
and the review fixes.

Removed: `transport_origin`/`to_origin`/`from_origin` and the origin
machinery, `NoTransport`, `transport_to_mvstd`, per-measure
`Base.rand(rng, T, μ)` methods, the `_trafo_cdf`/`_trafo_quantile` hooks.

Behavior changes for NEWS: univariate Distributions pivot on
`StdLogistic`; nested powers and array products return ArraysOfArrays
views, powers of tuple products struct arrays; `Half` transports via
`StdUniform` (tails limited); `mcombine(vcat, ::Product, ::Product)`
merges only concrete homogeneous marginals; `rand(rng, Int, μ)`
unsupported; isbits marginal arrays become struct arrays (`rand` of such
products gives plain arrays); vcat-combined and bind variates are flat;
`transport_to(StdUniform, m)` with binds inside errors (use
`StdUniform()^n`); batched kernels of rank-less measures throw.

## Verification

Full suite (Aqua, extensions, doctests) on CPU with JLArrays cases;
`test/reactant` (opt-in, backend via `MEASUREBASE_REACTANT_BACKEND`) and
`test/cuda` (opt-in) run locally on the GB10, both green at HEAD except
one expected-broken CUDA case (AffineMaps Jacobian on device).

## Known gaps and open decisions

- Upstream: AffineMaps lacks `Adapt` rules and device/traced Jacobians;
  Distributions isn't device-aware; ChangesOfVariables has no rules for
  `Base.Fix1`/`Fix2` arithmetic; HeterogeneousComputing has no Reactant
  compute unit; JLArrays has no RNG; Reactant rejects traced
  `VectorOfArrays` and empty batches.
- Decisions pending: one convention for out-of-support inputs (NaN mask
  vs. DomainError vs. AssertionError), `Half` tails via log-ccdf, device
  random variate infrastructure and `rand!`, Tier-1 static variates,
  the `smart-constructors.jl` review (location-scale arrays as affine
  pushforwards of powers), `_static_ndims` type-first vs. instance-first.
- Polish before merge: docs pass, NEWS, history curation, version bump,
  remove this file.
