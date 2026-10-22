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

@testset "Lebesgue measure on the standard simplex" begin
    using Random: Xoshiro
    using Statistics: mean, var
    using MeasureBase: StdSimplex, getdof
    @test [0.2, 0.3, 0.5] ∈ StdSimplex(3)
    @test [0.2, 0.3, 0.6] ∉ StdSimplex(3) && [-0.2, 0.7, 0.5] ∉ StdSimplex(3) && [0.5, 0.5] ∉ StdSimplex(3)
    L = Lebesgue(StdSimplex(3))
    @test basemeasure(L) === L && rootmeasure(L) === L
    @test getdof(L) == 2
    # The coordinate measure of the triangle with unit legs, 1/(n - 1)! in general:
    @test massof(L) == 1 / 2 && massof(Lebesgue(StdSimplex(5))) ≈ 1 / 24
    @test logdensityof(L, [0.2, 0.3, 0.5]) == 0 && logdensityof(L, [0.2, 0.3, 0.6]) == -Inf
    @test logdensities(L, [0.2 0.5 -0.1; 0.3 0.6 0.6; 0.5 0.1 0.5]) == [0, -Inf, -Inf]
    @test_throws ArgumentError logdensityof(L, [0.5, 0.5])
    # The normalized law is uniform on the simplex, Dirichlet(1, 1, 1):
    X = MeasureBase.batched_rand_impl(MeasureBase.GenContext{Float64}(Xoshiro(42)), L, (10^4,))
    @test size(X) == (3, 10^4) && all(x -> x ∈ StdSimplex(3), eachcol(X))
    @test vec(mean(X; dims = 2)) ≈ fill(1 / 3, 3) atol = 0.01
    @test vec(var(X; dims = 2)) ≈ fill(1 / 18, 3) atol = 0.005
    @test testvalue(L) ≈ fill(1 / 3, 3)
end
