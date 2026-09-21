# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Static variates end to end: measures whose variate sizes are statically
# known generate, transport and evaluate static arrays, type stable and
# allocation free.

using Test

using MeasureBase
using MeasureBase: StdNormal, StdUniform, StdExponential
using MeasureBase: productmeasure, mbind, mcombine, weightedmeasure, transport_to, logdensityof
using MeasureBase: batched_logdensityof_with_rest, batched_transport_to_std,
    batched_transport_from_std, batched_transport_to_std_with_rest,
    batched_transport_from_std_with_rest, _materialize
using MeasureBase.InverseFunctions: inverse
using ArraysOfArrays: flatview, sliced
using StaticArrays: SVector, SMatrix, Size
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

        @test allocations_of(rand, StdNormal()^static(3)) == 0

        # Nested powers keep the flat `(base dims..., power dims...)` rule:
        x = @inferred(rand((StdNormal()^static(2))^static(3)))
        @test flatview(x) isa SMatrix{2,3,Float64}
        @test length(x) == 3 && all(xi -> xi isa SVector{2,Float64}, x)
        @test reduce(hcat, x) == flatview(x)

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
        @test allocations_of(rand, μ) == 0

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
        @test allocations_of(rand, μ) == 0
        @test allocations_of(rand, productmeasure((a = StdNormal(), b = StdExponential()))) == 0

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
        ν = StdUniform()^static(3)
        X = SMatrix{3,4}(reshape(collect(1:12) ./ 10, 3, 4))

        @test @inferred(logdensities(μ, X)) ≈
              [logdensityof(μ, SVector{3}(X[:, i])) for i in 1:4]
        @test allocations_of(logdensities, μ, X) == 0

        Z = @inferred batched_transport_to_std(StdUniform, μ, X)
        @test Z ≈ reduce(hcat, [transport_to(ν, μ)(SVector{3}(X[:, i])) for i in 1:4])
        @test @inferred(batched_transport_from_std(StdUniform, μ, Z)) ≈ X

        # The broadcast hook transports the whole batch at once:
        Y = transport_to(ν, μ).(sliced(X, Val(1)))
        @test flatview(Y) ≈ Z
        @test Y[2] ≈ transport_to(ν, μ)(SVector{3}(X[:, 2]))
    end

    # Several variates per stream, with the multiplicity as a tuple of
    # static integers and as a `StaticArrays.Size`:
    @testset "static stream multiplicity" begin
        μ = StdNormal()^static(2)
        mc = mcombine(vcat, StdNormal()^static(2), StdUniform()^static(3))
        x = SVector{4}(randn(4))
        xc = SVector{10}(vcat(randn(2), rand(3), randn(2), rand(3)))
        to_u = transport_to(StdUniform(), StdNormal())

        for sz in ((static(2),), Size(2))
            z, x_rest = batched_transport_to_std_with_rest(StdUniform, μ, x, sz)
            @test z isa SVector{4,Float64} && isempty(x_rest)
            @test z ≈ to_u.(x)
            x_back, z_rest = batched_transport_from_std_with_rest(StdUniform, μ, z, sz)
            @test x_back ≈ reshape(x, (2, 2)) && size(z_rest, 1) == 0

            ℓ, x_ld_rest = batched_logdensityof_with_rest(μ, x, sz)
            @test _materialize(ℓ) isa SVector{2,Float64} && isempty(x_ld_rest)
            @test _materialize(ℓ) ≈ [logdensityof(μ, x[(2i - 1):(2i)]) for i in 1:2]

            zc, xc_rest = batched_transport_to_std_with_rest(StdUniform, mc, xc, sz)
            @test length(zc) == 10 && isempty(xc_rest)
            ℓc, xc_ld_rest = batched_logdensityof_with_rest(mc, xc, sz)
            @test isempty(xc_ld_rest)
            @test _materialize(ℓc) ≈ [logdensityof(mc, xc[(5i - 4):(5i)]) for i in 1:2]
        end
    end
end
