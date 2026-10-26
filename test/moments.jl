using Test
using LinearAlgebra: Diagonal, diag
using Statistics: mean, var, std, cov
using IntervalSets

using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential, StdLogistic, normalize
using AffineMaps: Mul, Add, MulAdd, AddMul, InvMul, Subtract, InvMulAdd, InvAddMul
using StaticArrays: SVector
using Static: static

@testset "moments" begin
    @testset "standard measures" begin
        for (μ, m, v) in [
            (StdUniform(), 1 / 2, 1 / 12),
            (StdNormal(), 0, 1),
            (StdExponential(), 1, 1),
            (StdLogistic(), 0, π^2 / 3),
        ]
            @test @inferred(mean(μ)) ≈ m
            @test @inferred(var(μ)) ≈ v
            @test @inferred(std(μ)) ≈ sqrt(v)
        end
    end

    @testset "point masses" begin
        d = MeasureBase.Dirac(1.5)
        @test mean(d) === 1.5 && var(d) === 0.0 && std(d) === 0.0
        @test var(MeasureBase.Dirac(2)) === 0.0
        dv = MeasureBase.Dirac([1, 2, 3])
        @test mean(dv) == [1, 2, 3]
        @test var(dv) == zeros(3) && eltype(var(dv)) == Float64
        @test cov(dv) == zeros(3, 3)
        @test var(MeasureBase.Dirac((a = 1.0, b = [2.0, 3.0]))) == (a = 0.0, b = [0.0, 0.0])
    end

    @testset "weighted and normalized measures" begin
        @test mean(2.5 * StdExponential()) == 1
        @test var(normalize(2.5 * StdUniform())) ≈ 1 / 12
        @test cov(0.3 * StdNormal()^3) == Diagonal(ones(3))
    end

    @testset "products and powers" begin
        μ = StdUniform()^(2, 3)
        @test mean(μ) == fill(0.5, 2, 3)
        @test var(μ) ≈ fill(1 / 12, 2, 3)
        @test std(μ) ≈ fill(sqrt(1 / 12), 2, 3)
        @test cov(μ) ≈ Diagonal(fill(1 / 12, 6))
        # Static sizes stay static:
        @test mean(StdNormal()^static(3)) isa SVector{3}

        p = productmeasure([StdNormal(), StdExponential(), 2 * StdUniform()])
        @test mean(p) ≈ [0, 1, 0.5]
        @test var(p) ≈ [1, 1, 1 / 12]
        @test cov(p) ≈ Diagonal([1, 1, 1 / 12])

        tp = productmeasure((a = StdExponential(), b = StdNormal()^2))
        @test mean(tp) == (a = 1.0, b = [0.0, 0.0])
        @test var(tp) == (a = 1.0, b = [1.0, 1.0])
        @test all(std(productmeasure((StdUniform(), StdLogistic()))) .≈ (sqrt(1 / 12), π / sqrt(3)))
        @test_throws ArgumentError cov(tp)

        # Powers of array-variate measures have nested moments:
        pp = (StdNormal()^2)^3
        @test mean(pp) == fill(zeros(2), 3)
        @test var(pp) == fill(ones(2), 3)
        @test_throws ArgumentError cov(pp)
    end

    @testset "affine pushforwards" begin
        # Scalar variates:
        for (f, m, v) in [
            (Mul(2.0), 1.0, 4 / 12),
            (Add(1.5), 2.0, 1 / 12),
            (Subtract(1.5), -1.0, 1 / 12),
            (MulAdd(-3.0, 1.0), -0.5, 9 / 12),
            (AddMul(1.0, 2.0), 3.0, 4 / 12),
            (InvMul(2.0), 0.25, 1 / 48),
            (InvMulAdd(2.0, 1.0), -0.25, 1 / 48),
            (InvAddMul(1.0, 2.0), -0.75, 1 / 48),
        ]
            ν = pushfwd(f, StdUniform())
            @test mean(ν) ≈ m
            @test var(ν) ≈ v
            @test std(ν) ≈ sqrt(v)
        end

        # Vector variates:
        A = [2.0 0.5; -1.0 1.5]
        b = [0.3, -0.7]
        μ = productmeasure([StdExponential(), 2 * StdUniform()])
        C = Diagonal([1, 1 / 12])
        for (f, m, Σ) in [
            (MulAdd(A, b), A * [1, 0.5] + b, A * C * A'),
            (AddMul(b, A), A * ([1, 0.5] + b), A * C * A'),
            (InvMul(A), A \ [1, 0.5], inv(A) * C * inv(A)'),
            (Add(b), [1, 0.5] + b, C),
        ]
            ν = pushfwd(f, μ)
            @test mean(ν) ≈ m
            @test cov(ν) ≈ Σ
            @test var(ν) ≈ diag(Σ)
        end
        @test mean(pushfwd(MulAdd(A, b), StdNormal()^2)) ≈ b
        @test cov(pushfwd(Mul(A), 3 * StdNormal()^2)) ≈ A * A'
    end

    @testset "unknown moments" begin
        for μ in [
            truncate(StdNormal(), 0 .. 1),
            Lebesgue(),
            restrict(x -> x > 0, StdNormal()),
            pushfwd(exp, StdNormal()),
            superpose(StdNormal(), StdUniform()),
            Half(StdNormal()),
        ]
            @test_throws ArgumentError mean(μ)
            @test_throws ArgumentError var(μ)
            @test_throws ArgumentError std(μ)
        end
        @test_throws ArgumentError cov(StdNormal())
        @test_throws ArgumentError cov(pushfwd(MulAdd(2.0, 1.0), StdNormal()))
        @test_throws ArgumentError cov(pushfwd(Mul(2.0), StdNormal()^(2, 2)))
    end
end
