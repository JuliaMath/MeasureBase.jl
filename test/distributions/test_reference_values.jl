# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Densities, derivatives and tails of the wrapped distribution families
# against independent references: analytic values and higher precision.

using Test

using LinearAlgebra
using Distributions
import AffineMaps
using ArraysOfArrays: nestedview, flatview
import ForwardDiff, Zygote, Mooncake
using SpecialFunctions: digamma

using MeasureBase
using MeasureBase: StdNormal, StdUniform, transport_to, logdensities, pushfwd

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

    @testset "shape derivatives of the incomplete gamma and beta functions" begin
        # Reference values from mpmath (numerical differentiation at 50 digits):
        for (a, x, ref) in [
            (0.01, 1.0e-5, -9.8167252611925644),
            (0.01, 0.5, -0.56554591171525781),
            (0.01, 20.0, -1.0564460222193207e-10),
            (0.5, 0.01, -0.52251433779748302),
            (0.5, 1.9, -0.15075191952357514),
            (2.5, 2.75, -0.245741695389706),
            (2.5, 125.0, -2.2719815141074742e-51),
            (100.0, 80.0, -0.0044026432434686354),
            (100.0, 117.0, -0.010053392530833614),
            (10000.0, 9800.0, -0.00053081913130245152),
            (10000.0, 10000.0, -0.0039894560494536683),
            (10000.0, 10300.0, -4.7684995140463611e-5),
        ]
            @test MeasureBase._gamma_cdf_dshape(a, x) ≈ ref rtol = 1e-10
        end
        for (a, b, x, ref_a, ref_b) in [
            (0.01, 0.01, 0.5, -25.004009427632699, 25.004009427632699),
            (0.01, 10.0, 1.0e-30, -34.160918876479137, 0.00054187794405317934),
            (0.1, 3.0, 0.2, -0.52003687114105701, 0.021625588481686938),
            (2.0, 3.0, 0.3, -0.2460937277467235, 0.12954506096487189),
            (2.0, 3.0, 0.7, -0.069571251131060496, 0.08187252372208087),
            (3.0, 0.1, 0.99, -0.027525027681462349, 2.3542697095364313),
            (30.0, 30.0, 0.45, -0.040413295224917996, 0.036568626794300736),
            (0.5, 10000.0, 1.0e-5, -0.82292693041303622, 1.6143697022234938e-5),
            (10000.0, 0.5, 0.99999, -1.6143697022205549e-5, 0.82292693041278708),
            (1000.0, 10000.0, 0.09, -0.011439648473102025, 0.001137496159575596),
            (10000.0, 10000.0, 0.497, -0.0019740649137160666, 0.0019622559785972475),
            (0.01, 10000.0, 0.0001, -0.22385735472191154, 3.699968863145381e-7),
        ]
            ∂a, ∂b = MeasureBase._beta_cdf_dshape(a, b, x)
            @test ∂a ≈ ref_a rtol = 1e-10
            @test ∂b ≈ ref_b rtol = 1e-10
        end
        @test MeasureBase._gamma_cdf_dshape(2.5f0, 2.75f0) ≈ -0.245741695389706 rtol = eps(Float32)
        ∂32 = MeasureBase._beta_cdf_dshape(2f0, 3f0, 0.3f0)
        @test ∂32 isa NTuple{2,Float32} && all(isapprox.(∂32, MeasureBase._beta_cdf_dshape(2.0, 3.0, Float64(0.3f0)), rtol = eps(Float32)))
        # Constant at the ends of the support, NaN outside of the domain:
        @test MeasureBase._gamma_cdf_dshape(2.0, 0.0) == 0 == MeasureBase._gamma_cdf_dshape(2.0, Inf)
        @test isnan(MeasureBase._gamma_cdf_dshape(2.0, -1.0)) && isnan(MeasureBase._gamma_cdf_dshape(0.0, 1.0))
        @test MeasureBase._beta_cdf_dshape(2.0, 3.0, 0.0) == (0, 0) == MeasureBase._beta_cdf_dshape(2.0, 3.0, 1.0)
        @test all(isnan, MeasureBase._beta_cdf_dshape(2.0, 3.0, 1.5)) && all(isnan, MeasureBase._beta_cdf_dshape(2.0, -1.0, 0.5))
    end

    @testset "shape derivatives of transports" begin
        # Fourth-order central differences as references:
        fd(f, a::Real) = (h = 1e-4 * a; (8 * (f(a + h) - f(a - h)) - (f(a + 2h) - f(a - 2h))) / (12 * h))
        fd(f, a::AbstractVector) = [fd(t -> f(setindex!(copy(a), t, i)), a[i]) for i in eachindex(a)]
        # With ForwardDiffPullbacks loaded, Zygote fails to pass the thunked
        # distribution tangents of its univariate transport rules to the
        # distribution constructors (for all distributions):
        fdpullbacks = !isnothing(Base.get_extension(MeasureBase, :MeasureBaseDistributionsForwardDiffPullbacksExt))
        @testset "$backend" for (backend, ∂) in pairs(_derivatives)
            for shape in [0.05, 0.8, 3.0, 50.0, 2000.0]
                x_gamma = quantile(Gamma(shape, 1.5), 0.3)
                x_beta, y_beta = quantile(Beta(shape, 2.0), 0.6), quantile(Beta(1.5, shape), 0.6)
                for f in [
                    a -> transport_to(StdUniform(), Gamma(a, 1.5))(x_gamma),
                    a -> transport_to(Gamma(a, 1.5), StdUniform())(0.97),
                    a -> transport_to(StdUniform(), Beta(a, 2.0))(x_beta),
                    b -> transport_to(StdUniform(), Beta(1.5, b))(y_beta),
                    a -> transport_to(Beta(a, 2.0), StdUniform())(0.2),
                    b -> transport_to(Beta(1.5, b), StdUniform())(0.9),
                ]
                    @test ∂(f, shape) ≈ fd(f, shape) rtol = 1e-7 broken = backend == :Zygote && fdpullbacks
                end
            end
            # Dirichlet transports, also with a concentration that depends on
            # another parameter:
            x, u, U = [0.2, 0.3, 0.5], [0.3, 0.8], [0.3 0.1 0.9; 0.8 0.5 0.2]
            for α in [[0.3, 2.0, 40.0], [2.0, 3.0, 4.0]]
                to_uniform = α -> dot([0.4, -0.7], transport_to(StdUniform()^2, Dirichlet(α))(x))
                from_uniform = α -> dot([0.2, -0.5, 0.7], transport_to(Dirichlet(α), StdUniform()^2)(u))
                @test ∂(to_uniform, α) ≈ fd(to_uniform, α) rtol = 1e-7
                @test ∂(from_uniform, α) ≈ fd(from_uniform, α) rtol = 1e-7
                hierarchical = τ -> from_uniform(τ .* α)
                @test ∂(hierarchical, 1.5) ≈ fd(hierarchical, 1.5) rtol = 1e-7
                batched = α -> sum(flatview(transport_to(Dirichlet(α), StdUniform()^2).(nestedview(U))) .* [0.2, -0.5, 0.7])
                @test ∂(batched, α) ≈ fd(batched, α) rtol = 1e-7
            end
            # Laws of transported variates:
            ℓ_gamma(a) = logdensityof(pushfwd(transport_to(Gamma(a, 1.5), StdNormal()), StdNormal()), 1.7)
            @test ∂(ℓ_gamma, 2.5) ≈ log(1.7 / 1.5) - digamma(2.5)
            α = [0.7, 2.0, 3.0]
            ℓ_dirichlet(α) = logdensityof(pushfwd(transport_to(Dirichlet(α), StdNormal()^2), StdNormal()^2), x)
            @test ∂(ℓ_dirichlet, α) ≈ log.(x) .- digamma.(α) .+ digamma(sum(α))
        end
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

        # Quantiles of the other family kernels (Distributions' Float64
        # quantiles serve as references, but not for Cauchy):
        ps = Float32[1f-37, 1f-20, 1f-9, 1f-3, 0.3, 0.5, 0.7, 1 - 1f-3, 1 - 6f-8]
        for (D, θ) in [(Normal, (0.3, 1.7)), (Exponential, (0.7,)), (Logistic, (0.2, 1.3)), (Laplace, (-0.4, 1.1)),
                       (LogNormal, (0.2, 0.6)), (Weibull, (1.4, 0.9)), (Gamma, (2.3, 1.2)), (Beta, (0.4, 0.6))]
            xs = transport_to(asmeasure(D(Float32.(θ)...)), StdUniform()).(ps)
            @test xs isa Vector{Float32} && issorted(xs)
            @test all(isapprox.(xs, Float32.(quantile.(D(θ...), Float64.(ps))), rtol = 16 * eps(Float32)))
        end

        # Logistic densities in both tails:
        d = Logistic(0.5f0, 2f0)
        for x in Float32[-floatmax(Float32), -1f30, -80, 0.25, 80, 1f30, floatmax(Float32)]
            z = abs(big(x) - d.μ) / d.θ
            @test logdensityof(asmeasure(d), x) ≈ Float32(-z - 2 * log1p(exp(-z)) - log(big(d.θ))) rtol = 4 * eps(Float32)
        end
    end
end
