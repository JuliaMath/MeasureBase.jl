# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Densities, derivatives and tails of the wrapped distribution families
# against independent references: analytic values and higher precision.

using Test

using LinearAlgebra
using Distributions
import AffineMaps
import ForwardDiff, Zygote, Mooncake
using SpecialFunctions: digamma

using MeasureBase
using MeasureBase: StdNormal, StdUniform, transport_to, logdensities

const _derivatives = (
    ForwardDiff = (f, x) -> x isa Real ? ForwardDiff.derivative(f, x) : ForwardDiff.gradient(f, x),
    Zygote = (f, x) -> something(only(Zygote.gradient(f, x)), zero(x)),
    Mooncake = (f, x) -> Mooncake.value_and_gradient!!(Mooncake.prepare_gradient_cache(f, x), f, x)[2][2],
)

@testset "reference values" begin
    @testset "densities of affine pushforwards" begin
        # The law of `2z` for a standard normal `z`:
        p = pushfwd(AffineMaps.Mul(2.0I(2)), StdNormal()^2)
        Y = [0.0 0.3 -1.2; 0.0 1.1 0.4]
        ℓ_ref = logpdf(MvNormal(zeros(2), 4.0I(2)), Y)
        @test logdensity_def(p, Y[:, 2]) ≈ ℓ_ref[2]
        @test MeasureBase.batched_logdensity_def(p, zeros(2)) ≈ ℓ_ref[1]
        @test collect(MeasureBase.batched_logdensity_def(p, Y)) ≈ ℓ_ref
        @test logdensities(p, Y) ≈ ℓ_ref
        @test logdensity_def(p^3, Y) ≈ sum(ℓ_ref)
        @test logdensity_rel(p^3, basemeasure(p)^3, Y) ≈ sum(ℓ_ref)
        @test logdensityof(p^3, Y) ≈ sum(ℓ_ref)
    end

    @testset "derivatives at unit shapes and medians" begin
        ld(d, x) = logdensityof(asmeasure(d), x)
        α, x = [1.0, 2.0, 3.0], [0.2, 0.3, 0.5]
        cases = [
            (a -> ld(Gamma(a, 1.5), 2.0), 1.0, log(2.0 / 1.5) - digamma(1.0)),
            (a -> ld(Beta(a, 2.0), 0.2), 1.0, log(0.2) - digamma(1.0) + digamma(3.0)),
            (b -> ld(Beta(2.0, b), 0.2), 1.0, log(0.8) - digamma(1.0) + digamma(3.0)),
            (a -> ld(Weibull(a, 1.5), 2.0), 1.0, 1 + log(2.0 / 1.5) * (1 - 2.0 / 1.5)),
            (a -> ld(Dirichlet(a), x), α, log.(x) .- digamma.(α) .+ digamma(sum(α))),
            (transport_to(asmeasure(Laplace(0.3, 1.3)), StdUniform()), 0.5, 2 * 1.3),
            (transport_to(StdUniform(), asmeasure(Laplace(0.3, 1.3))), 0.3, 1 / (2 * 1.3)),
            (transport_to(asmeasure(Cauchy(0.3, 1.3)), StdUniform()), 0.5, π * 1.3),
            (transport_to(StdUniform(), asmeasure(Cauchy(0.3, 1.3))), 0.3, 1 / (π * 1.3)),
            (transport_to(StdNormal(), asmeasure(Cauchy(0.3, 1.3))), 0.3, sqrt(2π) / (π * 1.3)),
        ]
        @testset "$backend" for (backend, ∂) in pairs(_derivatives)
            for (f, x, ref) in cases
                @test ∂(f, x) ≈ ref
            end
        end
        # The endpoints at unit shapes:
        @test ld(Gamma(1.0, 2.0), 0.0) ≈ -log(2.0)
        @test ForwardDiff.derivative(x -> ld(Gamma(1.0, 2.0), x), 0.0) ≈ -0.5
        @test ld(Beta(1.0, 2.0), 0.0) ≈ log(2.0) && ld(Beta(2.0, 1.0), 1.0) ≈ log(2.0)
        @test ld(Poisson(0.0), 0) == 0
        # Far tails, where the unselected branches overflow:
        @test _derivatives.Zygote(transport_to(StdUniform(), asmeasure(Laplace())), -800.0) == 0
        @test _derivatives.Zygote(transport_to(StdUniform(), asmeasure(Laplace(0f0, 1f0))), -90f0) ≈ exp(-90f0) / 2
    end

    @testset "Float32 tails" begin
        # Cauchy quantiles are `-cot(πp)`:
        f = transport_to(asmeasure(Cauchy(0f0, 1f0)), StdUniform())
        tail = Float32[10.0^-k for k in 37:-1:1]
        ps = vcat(tail, 0.3f0, 0.5f0, 0.7f0, 1 .- reverse(filter(>=(eps(Float32)), tail)))
        xs = f.(ps)
        @test xs isa Vector{Float32}
        @test all(isapprox.(xs, Float32.(-cot.(big(π) .* ps)), rtol = 4 * eps(Float32)))
        @test issorted(xs) && sign.(xs) == sign.(ps .- 0.5f0)
        @test f(1f-9) ≈ -1 / (π * 1f-9) rtol = 4 * eps(Float32)
        @test transport_to(asmeasure(Cauchy()), StdUniform())(1e-300) ≈ -1 / (π * 1e-300)

        # Cauchy densities where `z^2` overflows:
        d = Cauchy(0.5f0, 2f0)
        for x in Float32[1f20, -1f25, floatmax(Float32), -floatmax(Float32), 3, 0.25]
            z = (big(x) - d.μ) / d.σ
            @test logdensityof(asmeasure(d), x) ≈ Float32(-log1p(z^2) - log(big(π) * d.σ)) rtol = 4 * eps(Float32)
        end
        @test logdensityof(asmeasure(Cauchy(0f0, 1f0)), 1f20) ≈ logpdf(Cauchy(0.0, 1.0), 1e20)
    end
end
