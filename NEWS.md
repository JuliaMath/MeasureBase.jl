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
  - `StandardDist` (`StandardUniform`, `StandardNormal`) of the
    Distributions extension and the conversion of standard measures to
    distributions (`Distribution(StdNormal())`), use the standard measures
    `StdUniform()`, `StdNormal()` and their powers (as BAT 6 does).
    `transport_to(Normal, μ)` and `transport_to(Uniform, μ)` stay
    shorthands for `transport_to(StdNormal, μ)` and
    `transport_to(StdUniform, μ)`.
- `pullback` is deprecated in favor of `pullbck`.
- `insupport` and `params` are renamed to `insupportof` and `paramsof`,
  they clashed with the functions exported by Distributions. The old names
  are deprecated, measure implementations must add methods to the new
  names.
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
  never throw for values in the variate space of a measure. Pushforwards
  are only defined on the image of their map, outside of it the inverse
  may throw (`logdensityof(pushfwd(exp, StdNormal()), -1.0)` throws a
  `DomainError` from `log`).
- `getdof` may return a static integer, `insupportof` may return
  `MeasureBase.NoFastInsupport`.
- Batched kernels of measures without a declared variate rank throw.
- `==` and `hash` of measures are value-based.

### New features

- Batched densities: `logdensities(μ, X)` for flat batches, nested arrays,
  struct arrays and tuples of batches.
- Batched transports: `transports(f, X)` for `f = transport_to(ν, μ)`,
  with the batch conventions of `logdensities`. Broadcasting, `f.(X)`, is
  equivalent: it transports the batch as a whole, not variate by
  variate.
- Random variates via HeterogeneousComputing: `rand(ctx::GenContext, μ)`
  with RNG, precision and compute unit, batches via `rand(ctx, μ^n)`.
- Static variate sizes: statically sized measures give type stable,
  allocation-free densities, transports and draws.
- Arrays of isbits parameterized measures are stored as struct arrays and
  evaluated in fused kernels, also on GPUs and under Reactant.
- Parameter batches: broadcasting `logdensityof` and `densityof` over
  arrays of measures and over likelihoods (`logdensityof.(Ref(ℒ), Θ)` for
  arrays or struct arrays of parameter values) evaluates the measures in
  fused kernels, `logdensities(ms, X)` is the named form for flat batches.
  Density measures (posteriors) have batched kernels.
- `pushfwd(f, μ; varsize)` declares the variate size of a pushforward.
  `f` is no longer evaluated at construction, sizes that don't follow from
  types are learned from a test value when needed.
- `fast_dof` returns static or constant-time degrees of freedom and
  `MeasureBase.NoFastDOF` otherwise, transports and draws fall back to
  `getdof` or the with-rest path.
- Powers take integer-valued real counts (`StdNormal()^3.0`), so binds
  with counts drawn from discrete distributions work, e.g. the Poisson
  process `mbind(n -> μ^n, Poisson(λ), vcat)`.
- New exports: `logdensities`, `transports`, `mcombine`, `mbind`, `mkernel`,
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
  `MeasureBase.mspace_ndims`. The parts of binds, kernels and combined
  measures are available via `MeasureBase.kernelmeasure`,
  `MeasureBase.valuecombinator` and `MeasureBase.combinedmeasures`.
- Package extensions for Distributions (branch-free kernels for the main
  univariate families, `MvNormal`, `Dirichlet`), AffineMaps and Reactant.
- `truncate(μ, a..b)` restricts measures with scalar variates and no
  atoms to an interval (a `TruncatedMeasure`), structurally via their
  Stieltjes measure function `smf`: the mass is `smf(μ, b) - smf(μ, a)`,
  `normalize(truncate(μ, a..b))` is the conditional law, transports and
  draws go through `smf` and `invsmf` (tails limited by the resolution of
  uniform variates, as for `Half`). `smf` and `invsmf` are now defined for
  wrapped univariate distributions (cdf and quantile), weighted measures
  and `StdExponential`.
- `mean`, `var`, `std` and `cov` (Statistics) of the normalized laws of
  measures with exactly known moments: wrapped distributions, standard
  measures, point masses, weighted and normalized measures, products and
  powers (covariances for scalar marginals) and affine pushforwards
  (AffineMaps). Other measures throw an `ArgumentError`.
- The transports of `Gamma`, `Beta` and `Dirichlet` are differentiable
  with respect to their shape parameters under ForwardDiff, ChainRules
  (Zygote) and Mooncake, e.g. for hierarchical models with a Dirichlet
  concentration that depends on other parameters.

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
- `paramsof(::Type{<:ParameterizedMeasure})` (which returned the parameter
  names) is removed, use `paramnames`. Transformed measures have no
  parameters.
- Transports between point masses return floats.
- Numeric vectors and batches of tuple and named tuple products are flat
  streams, for densities and transports, like those of powers of tuple
  products (`[a, b1, b2]` for `(a, [b1, b2])`). Numeric vectors used to be
  read as vectors of marginal variates only, which is the same for scalar
  marginals.
- `testvalue` of wrapped univariate distributions is their median (the
  inversion of the constant uniform draw), also for discrete
  distributions like `DiscreteUniform` that failed before.
- Batches of point masses with numeric atoms (from `rand` and from
  transports, also as components of products and powers) are
  `FillArrays.Fill`s on the host, dense arrays on other compute units.
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
  Jacobians are density ratios of the normalized measures, they throw an
  `ArgumentError` where that requires an unknown mass.
- `normalize(μ)` returns a `NormalizedMeasure` (or `μ` itself if its mass
  is statically one), with a density if the mass of `μ` is known, it used
  to return a weighted measure. Normalized measures of unknown mass sample
  and transport, but have no density. `normalize` strips weights and throws
  for zero and infinite masses. `massof` is known for `Lebesgue(a..b)`,
  `Counting` of bounded integers, ranges and sets, and type stable for
  powers.
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

Measures with transports should declare `massof` (`static(1.0)` for
normalized measures), transport Jacobians of measures of unknown mass
throw.

### Dependencies

- New dependencies: HeterogeneousComputing (random variate generation),
  StaticThings, ArraysOfArrays, StructArrays, Adapt, ArgCheck and
  OneTwoMany.
- New weak dependencies: Distributions (with PDMats), AffineMaps,
  Reactant, ForwardDiff, ForwardDiffPullbacks and Mooncake.
- `Test` is a weak dependency instead of a dependency: the test utilities
  `MeasureBase.Interface.test_interface`, `test_transport` and `test_smf`
  are defined by a package extension and require `using Test`.
- Requires Compat 4.10 or later.
