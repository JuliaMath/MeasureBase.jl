# API

## Modules

```@index
Pages = ["api.md"]
Order = [:module]
```

## Types and constants

```@index
Pages = ["api.md"]
Order = [:type, :constant]
```

## Functions and macros

```@index
Pages = ["api.md"]
Order = [:macro, :function]
```

## Exported API

```@autodocs
Modules = [MeasureBase, MeasureBase.MeasureOperators]
Order = [:module, :type, :constant, :macro, :function]
Filter = t -> MeasureBaseDocs.api_section(t) == :exported
```

## Extension points

Public but not exported, extend and call them with the `MeasureBase.` prefix.

Measures with fixed-size variates declare `MeasureBase.mspace_ndims` (and
optionally `MeasureBase.mspace_flatsize`) and implement the batched
kernels. Measures with value-dependent variate sizes implement the
with-rest hooks (`MeasureBase.logdensityof_with_rest`,
`MeasureBase.transport_to_std_with_rest`,
`MeasureBase.transport_from_std_with_rest` and their `batched_` forms)
and declare `MeasureBase.fixed_stream_size`. They have no base measure,
so they also implement the point density and point transports,
`MeasureBase.rand_impl`, their degrees of freedom (`MeasureBase.NoDOF`)
and `MeasureBase.checked_arg`. See [Concepts](@ref) for both tiers, the
full lists and the stream conventions.

```@autodocs
Modules = [MeasureBase, MeasureBase.MeasureOperators]
Order = [:module, :type, :constant, :macro, :function]
Filter = t -> MeasureBaseDocs.api_section(t) == :public
```

Internal functions and types are documented under [Internals](@ref).

## Testing measure implementations

```@autodocs
Modules = [MeasureBase.Interface]
Order = [:module, :function]
```
