# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseMooncakeExt

using MeasureBase
import Mooncake
using Mooncake: @zero_derivative, MinimalCtx

using MeasureBase: isneginf, isposinf, _adignore_call
using MeasureBase: check_dof, require_insupport, _origin_depth
using MeasureBase: logdensityof_rt

# Unlike Zygote, Mooncake differentiates the collection utilities
# (`_pushfront`, etc., mutating code in general), `checked_arg` and
# `_checksupport` natively, so only the non-differentiable functions
# need rules:

@zero_derivative MinimalCtx Tuple{typeof(isneginf),Any}
@zero_derivative MinimalCtx Tuple{typeof(isposinf),Any}

@zero_derivative MinimalCtx Tuple{typeof(_adignore_call),Any}

@zero_derivative MinimalCtx Tuple{typeof(require_insupport),Any,Any}
@zero_derivative MinimalCtx Tuple{typeof(_origin_depth),Any}
@zero_derivative MinimalCtx Tuple{typeof(check_dof),Any,Any}

@zero_derivative MinimalCtx Tuple{typeof(logdensityof_rt),Any,Any}

end # module MeasureBaseMooncakeExt
