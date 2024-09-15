using DistributionMeasures
using Test

@testset "Distributions extension" begin
    include("test_autodiff_utils.jl")
    include("test_measure_interface.jl")
    include("test_distribution_measure.jl")
    include("test_standard_dist.jl")
    include("test_standard_uniform.jl")
    include("test_standard_normal.jl")
    include("test_transport.jl")
end
