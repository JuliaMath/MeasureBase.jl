# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseForwardDiffPullbacksExt

import MeasureBase
using ForwardDiffPullbacks: fwddiff

MeasureBase._fwddiff(f::Function) = fwddiff(f)

end # module MeasureBaseForwardDiffPullbacksExt
