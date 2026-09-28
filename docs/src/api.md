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

```@autodocs
Modules = [MeasureBase, MeasureBase.MeasureOperators]
Order = [:module, :type, :constant, :macro, :function]
Filter = t -> MeasureBaseDocs.api_section(t) == :public
```

Internal functions and types are documented under [Internals](@ref).
