# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Static variates end to end: measures whose variate sizes are statically
# known generate, transport and evaluate static arrays, type stable and
# allocation free.

using Test

using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential
using MeasureBase: productmeasure, mbind, weightedmeasure, transport_to, logdensityof
using MeasureBase.InverseFunctions: inverse
using ArraysOfArrays: ArrayOfSimilarArrays, flatview, sliced
using StaticArrays: SVector, SMatrix
using Static: static

include("testutils.jl")

# A hierarchical model over a named tuple of scalar and static-array
# marginals: the secondary marginals depend on the primary variate, their
# sizes don't.
const static_primary = productmeasure((
    a = StdNormal(),
    b = weightedmeasure(-0.5, StdExponential()),
))

static_kernel(x) = productmeasure((
    c = StdUniform()^static(2),
    d = weightedmeasure(-abs(x.a), StdNormal()^static(3)),
))

const static_model = mbind(static_kernel, static_primary, merge)

@testset "static variates" begin
    @testset "static powers" begin
        @test @inferred(rand(StdNormal()^static(3))) isa SVector{3,Float64}
        @test @inferred(rand(StdUniform()^static(2))) isa SVector{2,Float64}
        @test @inferred(rand(StdExponential()^static(4))) isa SVector{4,Float64}

        # Nested powers keep the flat `(base dims..., power dims...)` rule:
        x = @inferred(rand((StdNormal()^static(2))^static(3)))
        @test x isa ArrayOfSimilarArrays{Float64,1,1,<:SMatrix{2,3}}
        @test flatview(x) isa SMatrix{2,3,Float64}

        μ = StdNormal()^static(3)
        z = SVector(0.1, 0.2, 0.3)
        @test @inferred(transport_to(StdUniform()^static(3), μ)(z)) isa SVector{3,Float64}
        @test @inferred(transport_to(μ, StdUniform()^static(3))(SVector(0.1, 0.5, 0.9))) isa
              SVector{3,Float64}
        @test @inferred(logdensityof(μ, z)) ≈ sum(logdensityof.(Ref(StdNormal()), z))
        @test allocations_of(logdensityof, μ, z) == 0
        @test allocations_of(transport_to(StdUniform()^static(3), μ), z) == 0
    end

    @testset "hierarchical model with static marginals" begin
        μ = static_model
        x = @inferred rand(μ)
        @test x isa NamedTuple{(:a, :b, :c, :d)}
        @test x.a isa Float64
        @test x.b isa Float64
        @test x.c isa SVector{2,Float64}
        @test x.d isa SVector{3,Float64}

        ν = StdNormal()^static(7)
        f = transport_to(ν, μ)
        f_inv = inverse(f)

        z = @inferred f(x)
        @test z isa SVector{7,Float64}
        y = @inferred f_inv(z)
        @test y isa NamedTuple{(:a, :b, :c, :d)}
        @test y.c isa SVector{2,Float64}
        @test y.d isa SVector{3,Float64}
        @test all(map((u, v) -> u ≈ v, values(y), values(x)))

        @test allocations_of(f, x) == 0
        @test allocations_of(f_inv, z) == 0
        @test allocations_of(logdensityof, μ, x) == 0

        # The model density is the sum of the component densities:
        x_a = (a = x.a, b = x.b)
        x_b = (c = x.c, d = x.d)
        @test @inferred(logdensityof(μ, x)) ≈
              logdensityof(static_primary, x_a) + logdensityof(static_kernel(x_a), x_b)
    end

    @testset "products of scalar and static-array marginals" begin
        μ = productmeasure((
            a = StdNormal(),
            b = StdUniform()^static(2),
            c = weightedmeasure(-0.25, StdExponential()^static(3)),
        ))
        x = @inferred rand(μ)
        @test x isa NamedTuple{(:a, :b, :c)}
        @test x.a isa Float64
        @test x.b isa SVector{2,Float64}
        @test x.c isa SVector{3,Float64}

        f = transport_to(StdNormal()^static(6), μ)
        f_inv = inverse(f)
        z = @inferred f(x)
        @test z isa SVector{6,Float64}
        y = @inferred f_inv(z)
        @test y isa NamedTuple{(:a, :b, :c)}
        @test all(map((u, v) -> u ≈ v, values(y), values(x)))
        @test allocations_of(f, x) == 0
        @test allocations_of(f_inv, z) == 0
        @test allocations_of(logdensityof, μ, x) == 0

        # Tuple products behave the same way:
        μ_t = productmeasure((StdNormal(), StdUniform()^static(2)))
        x_t = @inferred rand(μ_t)
        @test x_t isa Tuple{Float64,SVector{2,Float64}}
        f_t = transport_to(StdNormal()^static(3), μ_t)
        @test @inferred(f_t(x_t)) isa SVector{3,Float64}
        @test @inferred(inverse(f_t)(f_t(x_t))) isa Tuple{Float64,SVector{2,Float64}}
    end

    @testset "batches of static variates" begin
        μ = StdNormal()^static(3)
        X = SMatrix{3,4}(reshape(collect(1:12) ./ 10, 3, 4))
        ℓ = logdensities(μ, X)
        @test ℓ ≈ [logdensityof(μ, SVector{3}(X[:, i])) for i in 1:4]

        ν = StdUniform()^static(3)
        Y = transport_to(ν, μ).(sliced(X, Val(1)))
        @test flatview(Y) ≈ reduce(hcat, [transport_to(ν, μ)(SVector{3}(X[:, i])) for i in 1:4])
    end
end
