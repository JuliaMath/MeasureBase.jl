# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsMooncakeExt

using MeasureBase
import Distributions
import Mooncake
using Mooncake: @zero_derivative, MinimalCtx

using Distributions: Distribution
using MeasureBase: _dist_params_numtype

# The distribution transports themselves need no rules here: Mooncake
# provides rules for Distributions and StatsFuns/SpecialFunctions, so it
# differentiates the cdf/quantile-based transports natively. The
# ForwardDiffPullbacks-based rules for `transport_def` are a
# Zygote/ChainRules pathway.

@zero_derivative MinimalCtx Tuple{typeof(_dist_params_numtype),Distribution}

end # module MeasureBaseDistributionsMooncakeExt
