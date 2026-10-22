# MeasureBase.jl release notes

## v0.15.0

A breaking release: MeasureBase is now batched first and runs on GPUs and
under Reactant. Densities, transports and random variates have one
implementation per measure type for single variates and batches. See the
"Concepts" page of the documentation for the design.

### Breaking changes

- **`⊙` now denotes the pullback** of a measure (`ν ⊙ f == pullbck(f, ν)`),
  it used to construct a pointwise product of measures
  (`PointwiseProductMeasure`, removed). Use `mintegrate(ℓ, prior)` for
  posteriors.
- Removed:
  - `transport_origin`, `to_origin`, `from_origin`, `NoTransport` and
    `transport_to_mvstd`, transport now goes through standard measures
    (see below).
  - `kleisli`, `AbstractTransitionKernel` and `↣`, kernels are plain
    functions (`mkernel`, `mbind`).
  - `rebase` and `commonbase`.
  - `likelihood_ratio` and `log_likelihood_ratio`.
  - Power-weighted measures (`↑`).
  - `ZeroSet`, `CodimOne`, `Simplex`, `ℝ₊` and `𝕀`.
  - Per-measure `rand(rng, T, μ)` methods and `rand(rng, Int, μ)`, see
    "Random variates" below. `rand!` for product measures.
  - The `_trafo_cdf`/`_trafo_quantile` hooks of the Distributions
    integration.
- `pullback` is deprecated in favor of `pullbck`.
- `mspace_elsize` is not exported anymore (it is still available as
  `MeasureBase.mspace_elsize`), `MeasureBase.mspace_flatsize` and
  `MeasureBase.mspace_ndims` replace it for kernels.
- Measure implementations extend new functions, see "Migration for measure
  implementations" below.

### Conventions

- Densities are `-Inf` outside of the support, including non-integer
  values of counting-based measures and infinite values.
- Transports return `NaN` outside of the support of the source measure,
  finite inputs give finite results (uniform inputs are clamped into the
  open unit interval).
- Only variates of the wrong shape throw (an `ArgumentError`), kernels
  never throw for values.
- `getdof` may return a static integer, `insupport` may return
  `MeasureBase.NoFastInsupport`.
- Batched kernels of measures without a declared variate rank throw.
- `==` and `hash` of measures are value-based.

### New features

- Batched densities: `logdensities(μ, X)` for flat batches, nested arrays,
  struct arrays and tuples of batches.
- Batched transports: `f.(X)` for `f = transport_to(ν, μ)`.
- Random variates via HeterogeneousComputing: `rand(ctx::GenContext, μ)`
  with RNG, precision and compute unit, batches via `rand(ctx, μ^n)`.
- Static variate sizes: statically sized measures give type stable,
  allocation-free densities, transports and draws.
- Arrays of isbits parameterized measures are stored as struct arrays and
  evaluated in fused kernels, also on GPUs and under Reactant.
- New exports: `logdensities`, `mcombine`, `mbind`, `mkernel`,
  `bindkernel`, `boundmeasure`, `mintegrate`, `mintegrate_exp`,
  `localmeasure`, `fast_dof`, `asmeasure`, `powermeasure`,
  `weightedmeasure`, `superpose`, `half`, `as_likelihood`,
  `likelihood_kernel`, `likelihood_obs`, `mdomain`, `valdomain`,
  `combinesets`, `⋄` (pushforward) and `▷` (bind). `⊗`, `∫`, `∫exp`, `𝒹`
  and `log𝒹` are still exported, all operators are defined in the
  submodule `MeasureBase.MeasureOperators`.
- Extension points and internals for measure implementations are declared
  `public` (on Julia >= 1.11), e.g. `MeasureBase.logdensityof_impl`,
  `MeasureBase.batched_rand_impl`, `MeasureBase.transport_to_std` and
  `MeasureBase.mspace_ndims`.
- Package extensions for Distributions (branch-free kernels for the main
  univariate families, `MvNormal`, `Dirichlet`), AffineMaps and Reactant.

### Behavior changes

- Univariate distributions without a dedicated kernel transport via
  `StdLogistic`, transports give numerically different (more
  tail-accurate) results.
- Nested powers and array products return ArraysOfArrays views, powers
  of tuple products return struct arrays. `rand` of products over isbits
  marginals gives plain arrays.
- Variates of `mcombine(vcat, ...)` and of binds are flat vectors.
- `mcombine(vcat, ::ProductMeasure, ::ProductMeasure)` merges only
  products of concrete, homogeneous marginals.
- `Half` transports via `StdUniform`, with limited tail accuracy. `Half`
  requires scalar variates and a parent that is symmetric about zero,
  known asymmetric parents throw at construction.
- The base measure of a `Dirac` is the point mass itself and the base of
  `Counting(S)` keeps the set `S`, so densities of superpositions of
  discrete and continuous components are taken relative to Lebesgue
  measure plus counting measure on the atoms.
- `mcombine(merge, ...)` throws on overlapping names.
- Transports between point masses return floats.
- Variates of the wrong shape throw `ArgumentError` everywhere, also for
  reshaped measures (`DimensionMismatch` before).
- `transport_to(StdUniform, μ)` for measures `μ` with binds inside
  throws, use `transport_to(StdUniform()^n, μ)`.
- Pushforwards of measures with counting root measures (point masses,
  counting measures, discrete distributions) push the root measure forward
  instead of applying a volume element, injective maps keep the masses of
  atoms (`logdensityof(pushfwd(Mul(2), Poisson(2)), 4) ==
  logdensityof(Poisson(2), 2)`).
- Dirichlet measures have the coordinate Lebesgue measure on the standard
  simplex, `Lebesgue(MeasureBase.StdSimplex(n))`, as base and root measure
  instead of Lebesgue measure on `ℝⁿ`, relative densities with respect to
  measures with densities on `ℝⁿ` throw.
- Dirichlet transports use FlatPPL's stick-breaking orientation, with the
  break fractions `z_i ~ Beta(α_i, Σ_{j>i} α_j)` and
  `x_i = z_i * (1 - Σ_{j<i} x_j)`. The standard uniform coordinates of
  Dirichlet variates differ from BAT 5.x and MeasureBase 0.14, they are
  complemented per coordinate (`u ↦ 1 - u`). Small break fractions are
  resolved exactly, components after a dominant one only to the resolution
  of uniform variates close to one.
- Transport Jacobians exclude weights: weighted measures transport like
  their base, map and log-abs-det Jacobian alike
  (`with_logabsdet_jacobian(transport_to(StdNormal(), 2 * StdNormal()), x) == (x, 0)`).
- The projections `mcombine(firstarg/secondarg, α, β)` and
  `mbind(f, α, firstarg)` keep the mass of the dropped measure as a weight.
  `rand` of a bind samples ancestrally from the normalized laws of its
  components.

### Migration for measure implementations

Implement

- `MeasureBase.logdensityof_impl` (or `logdensity_def` and
  `basemeasure`) instead of `logdensityof`, plus
  `MeasureBase.batched_logdensityof_impl` and `MeasureBase.mspace_ndims`
  for array-valued variates,
- `MeasureBase.preferred_stdmeasure`, `MeasureBase.transport_to_std` and
  `MeasureBase.transport_from_std` (plus their `batched_` forms for
  array-valued variates) instead of `transport_origin`, `to_origin` and
  `from_origin`,
- `MeasureBase.batched_rand_impl` (and optionally `MeasureBase.rand_impl`)
  instead of `rand`. Measures with a standard transport get draws via
  transport by default.

### Dependencies

- New dependencies: HeterogeneousComputing (random variate generation),
  StaticThings, ArraysOfArrays, StructArrays, Adapt, ArgCheck and
  OneTwoMany.
- New weak dependencies: Distributions (with StatsBase, StatsFuns and
  PDMats), AffineMaps, Reactant, ForwardDiff, ForwardDiffPullbacks and
  Mooncake.
