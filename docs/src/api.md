# API

## Modules

```@index
Pages = ["api.md", "extensions.md"]
Order = [:module]
```

## Types and constants

```@index
Pages = ["api.md", "extensions.md"]
Order = [:type, :constant]
```

## Functions and macros

```@index
Pages = ["api.md", "extensions.md"]
Order = [:macro, :function]
```

## Exported API

```@autodocs
Modules = [MeasureBase, MeasureBase.MeasureOperators]
Order = [:module, :type, :constant, :macro, :function]
Filter = t -> MeasureBaseDocs.api_section(t) == :exported
```

The public extension points and the test utilities for measure
implementations are documented under [Extension API](@ref), internal
functions and types under [Internals](@ref).
