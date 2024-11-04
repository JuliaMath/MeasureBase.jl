using Test

using MeasureBase
using MeasureBase: insupport as measure_insupport

using DensityInterface: logdensityof

@testset "primitive" begin
    for (m, x) in [
        (MeasureBase.LebesgueBase(), -1.0f0),
        (MeasureBase.LebesgueBase(), 1.0f0),
        (Lebesgue(), 1.0f0),
        (Lebesgue(), -1.0f0),
        (MeasureBase.CountingBase(), -1.0f0),
        (MeasureBase.CountingBase(), 1.0f0),
        (Counting(), 2),
        (Counting(), 2.0f0),
        (Counting(), 1.5f0),
        (Dirac(4.2), 4.2f0),
        (Dirac(4.2), -1.0f0),
        (Dirac([1, 2, 3]), [1, 2, 3]),
        (Dirac([4, 5]), [4, 5]),
    ]
        @testset "$(nameof(typeof(m)))" begin
            for x in [-Inf, -1.2, -1, 0, 0.0, 1 // 2, 0.5, 1, 1.0, 2, 2.3]
                if measure_insupport(m, x)
                    @test @inferred(logdensityof(m, x)) ≈ 0
                    if x isa Real
                        ld = logdensityof(m, x)
                        ld isa float(typeof(x))
                    end
                    @test MeasureBase.unsafe_logdensityof(m, x) ≈ 0
                else
                    @test @inferred(logdensityof(m, x)) ≈ -Inf
                end
            end
        end
    end
end
