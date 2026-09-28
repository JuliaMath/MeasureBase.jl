# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsMooncakeExt

using MeasureBase
import Distributions
import Mooncake
using Mooncake: @zero_derivative, @from_rrule, MinimalCtx
using Base: IEEEFloat

using Distributions: Distribution
using MeasureBase: _dist_params_numtype, _gamma_cdf, _gamma_quantile, _beta_cdf, _beta_quantile

# Mooncake differentiates the generic cdf/quantile-based transports
# natively. The incomplete gamma and beta functions and their inverses use
# the ChainRules rules of the Distributions/ChainRulesCore extension, so
# that all backends agree (shape parameter derivatives are NaN):

@zero_derivative MinimalCtx Tuple{typeof(_dist_params_numtype),Distribution}

@from_rrule MinimalCtx Tuple{typeof(_gamma_cdf),IEEEFloat,IEEEFloat}
@from_rrule MinimalCtx Tuple{typeof(_gamma_quantile),IEEEFloat,IEEEFloat}
@from_rrule MinimalCtx Tuple{typeof(_beta_cdf),IEEEFloat,IEEEFloat,IEEEFloat}
@from_rrule MinimalCtx Tuple{typeof(_beta_quantile),IEEEFloat,IEEEFloat,IEEEFloat}

end # module MeasureBaseDistributionsMooncakeExt
