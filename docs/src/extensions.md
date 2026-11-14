# Extension API

These names are public but not exported. Extend and call them with the
`MeasureBase.` prefix.

Measures with fixed-size variates declare `MeasureBase.mspace_ndims` (and
optionally `MeasureBase.mspace_flatsize`) and implement the batched
kernels. Measures with value-dependent variate sizes implement the
with-rest hooks (`MeasureBase.logdensityof_with_rest`,
`MeasureBase.transport_to_std_with_rest`,
`MeasureBase.transport_from_std_with_rest` and their `batched_` forms)
and declare `MeasureBase.fixed_stream_size`. They have no base measure,
their densities, point transports and degrees of freedom follow from the
with-rest hooks, and they implement `MeasureBase.rand_impl`. See
[Concepts](@ref) for both tiers, the full lists and the stream
conventions.

```@autodocs
Modules = [MeasureBase, MeasureBase.MeasureOperators]
Order = [:module, :type, :constant, :macro, :function]
Filter = t -> MeasureBaseDocs.api_section(t) == :public
```

## Testing measure implementations

```@autodocs
Modules = [MeasureBase.Interface]
Order = [:module, :function]
```
