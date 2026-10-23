# Concepts

MeasureBase provides one implementation per measure type for densities,
transports and random variates, for single variates and for batches, on
the CPU, on GPUs and under Reactant.

## Batched first

A single variate is a batch with zero batch dimensions. The point API
(`logdensityof`, `transport_to(ν, μ)(x)`, `rand`) and the batched API
(`logdensities`, `f.(X)`, `rand(μ^n)`) run the same kernels, so they can't
drift apart. Static arrays keep the point path allocation-free.

```julia
using MeasureBase

μ = StdNormal()^3
X = rand(μ^100)            # batch of 100 variates, nested view of a 3×100 matrix
logdensities(μ, X)         # 100 log-densities
f = transport_to(StdUniform()^3, μ)
f.(X)                      # the batch transported to the unit cube
```

**Variate rank.** [`MeasureBase.mspace_ndims`](@ref) gives the number of
dimensions of the variates of a measure type, `0` for scalar variates. It
is derived from [`MeasureBase.mspace_flatsize`](@ref) where the size is
known, declared by array-variate leaf measures and derived by structural
measures. Kernels only need the rank. Sizes are optional declarations,
used for validation, for stream consumption and for static fast paths,
never for routing. Without a rank, the batched defaults throw with a
message naming the missing declaration, point kernels keep working.

**Flat batches.** A flat batch is an array `(variate dims..., batch
dims...)`. Batched kernels return arrays over the batch dimensions (a
number for a single variate), possibly lazily. Entry points accept nested
arrays (ArraysOfArrays containers are fused into their flat storage,
ragged containers are evaluated variate by variate), tuples and named
tuples of batches for tuple-valued variates, struct arrays and arrays of
tuples. `rand` of a power of a tuple product returns a struct array.

**Static variates.** Sizes are `StaticThings.SizeLike` throughout. On the
CPU, bulk draws of a fully static size are static arrays and stream splits
and reshapes keep them static: `rand`, `transport_to` and `logdensityof`
of a measure whose variate sizes are all static are type stable and
allocation-free, for powers, tuple and named tuple products and binds
over them. Other compute units allocate their own arrays.

## Support and shape conventions

- Densities are `-Inf` outside of the support, including non-integer
  values of measures over counting measures and infinite values. `NaN`
  inputs give `NaN` or `-Inf`. Relative densities keep `+Inf`, `-Inf` and
  `NaN`.
- Transports return `NaN` outside of the support of the source measure.
  Finite inputs give finite results: uniform inputs are clamped into the
  open unit interval before quantiles, tail probabilities of log-space
  conversions floor at the smallest positive float.
- Variates of the wrong shape (rank, size, container or element kind,
  tuple arity) throw an `ArgumentError` at the entry points, raised by
  [`MeasureBase.checked_arg`](@ref) and, for densities, by comparing flat
  batches with [`MeasureBase.mspace_flatsize`](@ref) where it is known.
  Variates of the right shape never throw.
- Kernels are branch-free where a condition may be a device or traced
  value (support masks evaluate both branches), so implementations must
  not throw outside the support.

## Streams

Variates of `mcombine(vcat, ...)`, of binds and of tuple products inside
such combinations are flat vectors, consumed with the with-rest protocol.
Point forms take a single stream `x` (a vector) and return
`(result, x_μ, x_rest)`, with the consumed variate `x_μ` (binds need it)
and the unconsumed rest of the stream.

Batched forms take a batch of streams `X`, an array `(rows, batch
dims...)` with the streams along the first dimension, and a multiplicity
`sz::SizeLike`: each stream holds `prod(sz)` consecutive variates of the
measure, which are consumed together, and they return `(result, X_rest)`.
The result is a batch over `(sz..., batch dims...)`: log-densities as an
array of that shape, variates as a flat batch `(variate dims..., sz...,
batch dims...)`, standard variates as streams `(dof * prod(sz), batch
dims...)`. Scalar variates are consumed as `(1, sz..., batch dims...)`.
A single stream with `sz == ()` is the point case. Powers pass their size
as multiplicity to their base, combined measures and tuple products split
rows by their fixed stream lengths.

[`MeasureBase.fixed_stream_size`](@ref) decides whether a batch of
streams is consumed in fused operations or stream by stream: measures for
which it is false are only called with single streams and `sz == ()`, the
outermost stream combinator loops over the streams.

## Transport

`transport_to(ν, μ)` returns a [`MeasureBase.TransportFunction`](@ref).
`f(x)` transports a variate, `f.(X)` a whole batch. Transports are
invertible (`InverseFunctions.inverse`) and provide their log-abs-det
Jacobian (`ChangesOfVariables.with_logabsdet_jacobian`).

Transport goes through standard measures (`StdUniform`, `StdNormal`,
`StdExponential`, `StdLogistic` and their powers) as pivots. Each measure
type states its preferred standard measure
([`MeasureBase.preferred_stdmeasure`](@ref)), the pivot for a pair of
measures is chosen by [`MeasureBase.promote_stdmeasure`](@ref), and
[`MeasureBase.stdconvert`](@ref) converts between standard measures in
log form, tail-accurately. Measure types implement transport to and from
their preferred standard measure only. Standard variates are flat
streams `(dof, batch dims...)`.

[`getdof`](@ref) gives the degrees of freedom of a measure, and
[`fast_dof`](@ref) gives them where they are static or take constant time:
for standard measures, Lebesgue measure, point masses, wrapped univariate
and Dirichlet distributions, and for powers, tuple products, products over
arrays of marginals with one degree of freedom each and combinations of
measures with fast degrees of freedom. Other measures return
`MeasureBase.NoFastDOF`. Transports and random variate generation split
streams of standard variates into the chunks of single variates with
`fast_dof`. Without it, composed measures consume their standard variates
component by component (the with-rest path, see [Streams](@ref)) and
other measures fall back to `getdof`, once per call. Neither is used
inside kernels.

## Random variates

`rand(ctx::GenContext, μ)` draws a variate with the random number
generator, precision and compute unit of the generative context `ctx`
(from HeterogeneousComputing), `rand(μ)`, `rand(rng, μ)` and `rand(T, μ)`
are shorthands. Batches of variates are drawn as variates of powers,
`rand(ctx, μ^n)`. By default, standard variates are drawn in bulk on the
compute unit and transported to `μ`.

Floating-point variates have the precision of the context, also through
transformations and combinators. Point masses generate their atom
exactly, and integer variates (e.g. of counting measures) stay integers.

## Measure combinators

- `powermeasure(μ, dims)` or `μ^dims`: variates `(base dims..., power
  dims...)`, the innermost base dimensions first. Results are nested views
  over the flat storage.
- `productmeasure(μs)` for arrays, tuples and named tuples of measures.
  Arrays of isbits parameterized marginals are stored as struct arrays and
  evaluated in fused kernels, also on GPUs and under Reactant.
- `mcombine(vcat, α, β)` and `mcombine(merge, α, β)`, `mbind(k, μ)`
  (flat stream variates for `vcat`, merged named tuples for `merge`).
- `pushfwd(f, μ)` and `pullbck(f, μ)`, with batched application for
  `Base.BroadcastFunction` and AffineMaps functions. The variate size of a
  pushforward is declared (`pushfwd(f, μ; varsize)`), the size of the
  origin for `Base.BroadcastFunction`s, or follows from the type `f`
  returns for variates of the variate rank (or static size) of the origin.
  Otherwise it is learned from `f` applied to a test value of `μ` when an
  operation needs it. `f` is never evaluated at construction.
- `weightedmeasure`, `superpose`, `restrict`, `mintegrate`,
  `mintegrate_exp`, `half`.

The operators `⋄` (pushforward), `⊙` (pullback), `▷` (bind), `⊗`
(product), `∫`, `∫exp`, `𝒹` and `log𝒹` are shorthands for these, see
[`MeasureBase.MeasureOperators`](@ref).

## Extension points

Measure types implement the following, all other functionality is
derived.

| Aspect | Scalar-variate leaf | Array-variate leaf | Composed measure |
|:--- |:--- |:--- |:--- |
| Density | `logdensity_def` (or `logdensityof_impl`), `basemeasure`, `insupportof` | `mspace_ndims`, `batched_logdensityof_impl` (or `batched_logdensity_def`) | kernels via components, with-rest forms |
| Transport | `preferred_stdmeasure`, `transport_to_std`, `transport_from_std` | also `batched_transport_to_std`, `batched_transport_from_std` | with-rest forms, `fixed_stream_size` |
| Random | `batched_rand_impl` or `rand_impl` (derived from a transport) | `batched_rand_impl` | derived |
| Declarations | `mspace_flatsize` | `mspace_ndims`, optionally `mspace_flatsize` | derived |
| Degrees of freedom | `getdof` and `fast_dof` where they differ from those of the base measure (e.g. over counting measures) | the same | derived |

All of these are declared `public` (see the [API](@ref) page) and
must be qualified with `MeasureBase.` when extended.

Measures come in two tiers:

- **Fixed-size variates.** Measures whose variates all have the same size
  declare their variate rank, `MeasureBase.mspace_ndims(::Type{<:M})`, and
  optionally their flat variate size, `MeasureBase.mspace_flatsize(::M)`
  (or the type form, for static sizes). They get the fused batched path:
  batches of densities, transports and draws, and the consumption of any
  number of variates per stream in a batch of streams, from their batched
  kernels. Their fast degrees of freedom (`fast_dof`) let transports and
  draws split batches of standard variates up front.
- **Value-dependent variate sizes.** Measures whose variate size depends on
  the variate (like binds) implement the with-rest hooks
  [`MeasureBase.logdensityof_with_rest`](@ref),
  [`MeasureBase.batched_logdensityof_with_rest`](@ref),
  [`MeasureBase.transport_to_std_with_rest`](@ref),
  [`MeasureBase.transport_from_std_with_rest`](@ref),
  [`MeasureBase.batched_transport_to_std_with_rest`](@ref) and
  [`MeasureBase.batched_transport_from_std_with_rest`](@ref) (the batched
  forms for single streams and `sz == ()`), and declare
  `MeasureBase.fixed_stream_size(::Type{<:M}) = static(false)`. Batches
  of such variates are evaluated variate by variate, see
  [Streams](@ref) for the conventions.

A scalar-variate leaf `M` with `logdensity_def`, `basemeasure` and
`insupportof` has point densities. Batched densities, powers and products
need its variate rank, `MeasureBase.mspace_flatsize(::Type{<:M}) = ()`.
Streams (`mcombine(vcat, ...)`, binds) need its variate size,
`MeasureBase.mspace_flatsize(::M) = ()`, or draw a test variate to learn
it. Random variates need `batched_rand_impl` or `rand_impl`, or the
standard transport (`preferred_stdmeasure`, `transport_to_std`,
`transport_from_std`), which also provides `transport_to`.

- [`MeasureBase.logdensityof_impl`](@ref) implements `logdensityof`,
  [`MeasureBase.batched_logdensityof_impl`](@ref) implements
  `logdensities`. Measure types should specialize these (or
  `logdensity_def` and `basemeasure`), not `logdensityof` itself.
- [`MeasureBase.batched_rand_impl`](@ref) is the primary extension point
  for random variates. [`MeasureBase.rand_impl`](@ref) defaults to it (with
  zero batch dimensions) and may be specialized for single variates. The
  default `batched_rand_impl` transports standard variates if the measure
  has a standard transport and calls `rand_impl` variate by variate
  otherwise, so measures that specialize `rand_impl` and have a standard
  transport should specialize `batched_rand_impl` as well. Structural
  measures draw single variates of their components via `rand_impl`.
- [`MeasureBase.tpmeasure_split_combined`](@ref) splits variates of
  binds and combinations with other value combination functions than
  `tuple`, `Pair`, `vcat` and `merge`.
- [`MeasureBase.transport_def`](@ref) may be specialized for direct
  transports between specific pairs of measure types.
