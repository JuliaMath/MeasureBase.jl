# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test
using MeasureBase
using Distributions
import ForwardDiff, ForwardDiffPullbacks, ChainRulesCore

const MeasureBaseDistributionsExt = Base.get_extension(MeasureBase, :MeasureBaseDistributionsExt)
@test MeasureBaseDistributionsExt isa Module

using .MeasureBaseDistributionsExt:
    StandardDist, StandardUniform, StandardNormal, DistributionMeasure, nonstddist

@testset "Distributions extension" begin
    include("test_autodiff_utils.jl")
    include("test_measure_interface.jl")
    include("test_distribution_measure.jl")
    include("test_standard_dist.jl")
    include("test_standard_uniform.jl")
    include("test_standard_normal.jl")
    include("test_conversions.jl")
    include("test_transport.jl")
    include("test_mooncake.jl")
end
