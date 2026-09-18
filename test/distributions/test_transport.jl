# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using Test

using LinearAlgebra
using InverseFunctions, ChangesOfVariables
using Distributions, ArraysOfArrays
using ArraysOfArrays: sliced, flatview
using StableRNGs
import ForwardDiff, Zygote
import PDMats

using MeasureBase: transport_to
using MeasureBase: StdUniform, StdNormal, StdExponential, StdLogistic
using .MeasureBaseDistributionsExt: _trafo_logcdf, _trafo_logccdf, _trafo_quantile, _trafo_cquantile

include("getjacobian.jl")


@testset "test_distribution_transform" begin
    function test_back_and_forth(trg, src)
        @testset "transform $(typeof(trg).name) <-> $(typeof(src).name)" begin
            x = rand(src)
            y = transport_to(trg, src)(x)
            src_v_reco = transport_to(src, trg)(y)

            @test x ≈ src_v_reco

            f = x -> transport_to(trg, src)(x)
            ref_ladj = logpdf(src, x) - logpdf(trg, y)
            @test ref_ladj ≈ logabsdet(getjacobian(f, x))[1]
        end
    end

    reshaped_rand(d::Distribution{Univariate}, n) = rand(d, n)
    reshaped_rand(d::Distribution{Multivariate}, n) = sliced(rand(d, n))

    function test_dist_trafo_moments(trg, src)
        unshaped(x) = first(torv_and_back(x))
        @testset "check moments of trafo $(typeof(trg).name) <- $(typeof(src).name)" begin
            X = reshaped_rand(src, 10^5)
            Y = transport_to(trg, src).(X)
            Y_ref = reshaped_rand(trg, 10^6)
            @test isapprox(mean(unshaped.(Y)), mean(unshaped.(Y_ref)), rtol = 0.5)
            @test isapprox(cov(unshaped.(Y)), cov(unshaped.(Y_ref)), rtol = 0.5)
        end
    end

    @testset "transforms-tests" begin
        stduvuni = StandardDist{Uniform}()
        stduvnorm = StandardDist{Normal}()

        uniform1 = Uniform(-5.0, -0.01)
        uniform2 = Uniform(0.01, 5.0)

        normal1 = Normal(-10, 1)
        normal2 = Normal(10, 5)

        stdmvnorm1 = StandardDist{Normal}(1)
        stdmvnorm2 = StandardDist{Normal}(2)

        stdmvuni2 = StandardDist{Uniform}(2)

        standnorm2_reshaped = reshape(stdmvnorm2, 1, 2)

        mvnorm = MvNormal([0.3, -2.9], [1.7 0.5; 0.5 2.3])
        beta = Beta(3,1)
        gamma = Gamma(0.1,0.7)
        dirich = Dirichlet([0.1,4])

        test_back_and_forth(stduvuni, stduvuni)
        test_back_and_forth(stduvnorm, stduvnorm)
        test_back_and_forth(stduvuni, stduvnorm)
        test_back_and_forth(stduvnorm, stduvuni)

        test_back_and_forth(stdmvuni2, stdmvuni2)
        test_back_and_forth(stdmvnorm2, stdmvnorm2)
        test_back_and_forth(stdmvuni2, stdmvnorm2)
        test_back_and_forth(stdmvnorm2, stdmvuni2)

        test_back_and_forth(beta, stduvnorm)
        test_back_and_forth(gamma, stduvnorm)
        test_back_and_forth(gamma, beta)

        test_back_and_forth(mvnorm, stdmvuni2)
        test_back_and_forth(stdmvuni2, mvnorm)

        test_back_and_forth(mvnorm, standnorm2_reshaped)
        test_back_and_forth(standnorm2_reshaped, mvnorm)
        test_back_and_forth(stdmvnorm2, standnorm2_reshaped)
        test_back_and_forth(standnorm2_reshaped, standnorm2_reshaped)

        test_dist_trafo_moments(normal2, normal1)
        test_dist_trafo_moments(uniform2, uniform1)

        test_dist_trafo_moments(beta, stduvnorm)
        test_dist_trafo_moments(gamma, stduvnorm)

        test_dist_trafo_moments(mvnorm, stdmvnorm2)
        test_dist_trafo_moments(dirich, stdmvnorm1)

        let
            mvuni = product_distribution([Uniform(), Uniform()])

            x = rand()
            @test_throws ArgumentError transport_to(stduvnorm, mvnorm)(x)
            @test_throws ArgumentError transport_to(stduvnorm, stdmvnorm1)(x)
            @test_throws ArgumentError transport_to(stduvnorm, stdmvnorm2)(x)

            x = rand(2)
            @test_throws ArgumentError transport_to(stduvnorm, mvnorm)(x)
            @test_throws ArgumentError transport_to(stduvnorm, stdmvnorm1)(x)
            @test_throws ArgumentError transport_to(stduvnorm, stdmvnorm2)(x)
        end
    end

    @testset "Custom cdf and quantile for dual numbers" begin
        Dual = ForwardDiff.Dual
        dual_normal = Normal(Dual(0, 1, 0, 0), Dual(1, 0, 1, 0))
        dual_x = Dual(0.5, 0, 0, 1)
        dual_p = Dual(0.3, 0, 0, 1)

        @test isapprox(_trafo_logcdf(dual_normal, dual_x), logcdf(dual_normal, dual_x), rtol = 10^-6)
        @test isapprox(_trafo_logcdf(Normal(0, 1), Dual(0.5, 1)), logcdf(Normal(0, 1), Dual(0.5, 1)), rtol = 10^-6)
        @test isapprox(_trafo_logccdf(dual_normal, dual_x), logccdf(dual_normal, dual_x), rtol = 10^-6)
        @test isapprox(_trafo_logccdf(Normal(0, 1), Dual(0.5, 1)), logccdf(Normal(0, 1), Dual(0.5, 1)), rtol = 10^-6)

        @test isapprox(_trafo_quantile(Normal(0, 1), Dual(0.3, 1)), quantile(Normal(0, 1), Dual(0.3, 1)), rtol = 10^-6)
        @test isapprox(_trafo_quantile(dual_normal, dual_p), quantile(dual_normal, dual_p), rtol = 10^-6)
        @test isapprox(_trafo_cquantile(Normal(0, 1), Dual(0.3, 1)), cquantile(Normal(0, 1), Dual(0.3, 1)), rtol = 10^-6)
        @test isapprox(_trafo_cquantile(dual_normal, dual_p), cquantile(dual_normal, dual_p), rtol = 10^-6)

        # Distributions whose cdf doesn't support dual numbers natively:
        beta = Beta(2.0, 3.0)
        dlogitcdf(d, x) = pdf(d, x) / (cdf(d, x) * ccdf(d, x))
        @test ForwardDiff.derivative(x -> transport_to(StdLogistic(), beta)(x), 0.3) ≈ dlogitcdf(beta, 0.3)
        x_b = transport_to(beta, StdLogistic())(-0.4)
        @test ForwardDiff.derivative(l -> transport_to(beta, StdLogistic())(l), -0.4) ≈ inv(dlogitcdf(beta, x_b))
    end

    @testset "tails of univariate transports" begin
        # Bounded and heavy-lower-tailed distributions lose the lower tail in
        # their quantile functions, so the ranges differ:
        for (d, ls) in [
            (Normal(0.3, 1.7), [-700.0, -40.0, -8.0, 0.0, 8.0, 40.0, 700.0]),
            (Weibull(0.7, 1.3), [-40.0, -8.0, 0.0, 8.0, 40.0, 700.0]),
            (truncated(Normal(0.2, 1.1), -3.0, 2.5), [-8.0, 0.0, 8.0]),
        ]
            for l in ls
                x = transport_to(d, StdLogistic())(l)
                @test insupport(d, x)
                @test isapprox(transport_to(StdLogistic(), d)(x), l, rtol = 1e-6, atol = 1e-12)
            end
        end
        for z in [-8.0, 8.0, 37.0]
            x = transport_to(Weibull(0.7, 1.3), StdNormal())(z)
            @test isfinite(x) && x > 0
            @test transport_to(StdNormal(), Weibull(0.7, 1.3))(x) ≈ z rtol = 1e-6
        end
    end

    @testset "trafo autodiff pullbacks" begin
        x = [0.6, 0.7, 0.8, 0.9]
        f = transport_to(Dirichlet([3.0, 4.0, 5.0, 6.0, 7.0]), Uniform)
        @test isapprox(ForwardDiff.jacobian(f, x), Zygote.jacobian(f, x)[1], rtol = 10^-4)
        f = inverse(transport_to(Normal, Dirichlet([3.0, 4.0, 5.0, 6.0, 7.0])))
        @test isapprox(ForwardDiff.jacobian(f, x), Zygote.jacobian(f, x)[1], rtol = 10^-4)
    end


    @testset "transport_to autosel" begin
        for (M,R) in [
            (StandardDist{Normal}, StandardDist{Normal})
            (Normal, StandardDist{Normal})
            (StandardDist{Uniform}, StandardDist{Uniform})
            (Uniform, StandardDist{Uniform})
        ]
            @test @inferred(transport_to(M, Weibull())) == transport_to(R(), Weibull())
            @test @inferred(transport_to(Weibull(), M)) == transport_to(Weibull(), R())
            @test @inferred(transport_to(M, MvNormal(float(I(5))))) == transport_to(R(5), MvNormal(float(I(5))))
            @test @inferred(transport_to(MvNormal(float(I(5))), M)) == transport_to(MvNormal(float(I(5))), R(5))
            @test @inferred(transport_to(M, StdExponential()^(2,3))) == transport_to(R(6), StdExponential()^(2,3))
            @test @inferred(transport_to(StdExponential()^(2,3), M)) == transport_to(StdExponential()^(2,3), R(6))
        end
    end

    @testset "affine transformed distributions" begin
        d = 2.0 * Weibull(0.7) + 1.0
        x = rand(StableRNG(789990641), d)
        u = transport_to(StdUniform(), d)(x)
        @test u ≈ cdf(d, x)
        @test transport_to(d, StdUniform())(u) ≈ x
        test_back_and_forth(StandardDist{Normal}(), d)
    end

    @testset "truncated distributions" begin
        d = truncated(Normal(0.3, 1.2), -0.5, 1.5)
        for u in [0.0, 0.25, 0.75, 1.0, prevfloat(1.0)]
            x = transport_to(d, StdUniform())(u)
            @test minimum(d) <= x <= maximum(d)
        end
        test_back_and_forth(StandardDist{Uniform}(), d)
    end

    @testset "products of distributions" begin
        pd = product_distribution([Weibull(0.7), Exponential(1.3), Normal(0.5, 2.0)])
        m = MeasureBase.asmeasure(pd)
        x = rand(StableRNG(789990641), pd)
        for trg in [StdUniform()^3, StdNormal()^3]
            y = transport_to(trg, m)(x)
            y_ref = map((d_i, x_i) -> transport_to(trg.parent, d_i)(x_i), pd.v, x)
            @test y ≈ y_ref
            @test transport_to(m, trg)(y) ≈ x
        end

        pd2 = product_distribution([Normal(2.0, 0.5), Weibull(1.2), Uniform(-1.0, 3.0)])
        m2 = MeasureBase.asmeasure(pd2)
        y = transport_to(m2, m)(x)
        @test transport_to(m, m2)(y) ≈ x
    end

    @testset "batched transport" begin
        mvn = MvNormal([0.3, -2.9], [1.7 0.5; 0.5 2.3])
        f = transport_to(StdNormal()^2, mvn)
        X = rand(StableRNG(789990641), mvn, 6)
        Y = f.(sliced(X, Val(1)))
        @test flatview(Y) ≈ stack(map(f, eachcol(X)))
        @test flatview(inverse(f).(Y)) ≈ X
        g = transport_to(StdNormal(), Weibull(0.7, 1.3))
        x = rand(StableRNG(789990641), Weibull(0.7, 1.3), 10)
        @test g.(x) ≈ map(g, x)
        pd = product_distribution([Weibull(0.7), Exponential(1.3), Normal(0.5, 2.0)])
        h = transport_to(StdNormal()^3, asmeasure(pd))
        Xp = rand(StableRNG(789990641), pd, 5)
        @test stack(h.(sliced(Xp, Val(1)))) ≈ stack(map(h, eachcol(Xp)))
        pn = product_distribution([Normal(1.0, 2.0), Normal(0.0, 3.0), Normal(2.0, 1.0)])
        mn = asmeasure(pn)
        @test MeasureBase.mspace_flatsize(mn) == (3,)
        hn = transport_to(StdUniform()^3, mn)
        Xn = rand(StableRNG(789990641), pn, 4)
        Yn = hn.(sliced(Xn, Val(1)))
        @test flatview(Yn) ≈ stack(map(hn, eachcol(Xn)))
        @test flatview(inverse(hn).(Yn)) ≈ Xn
        @test eltype(rand(StableRNG(1), Float32, mn)) == Float32
        @test eltype(flatview(rand(StableRNG(1), Float32, mn^3))) == Float32
        @test eltype(rand(StableRNG(1), Float32, asmeasure(pd))) == Float32
    end

    @testset "MvNormal covariance representations" begin
        for Σ in [PDMats.ScalMat(3, 2.5), PDMats.PDiagMat([0.5, 1.0, 2.5]), Diagonal([0.5, 1.0, 2.5])]
            mvn = MvNormal([0.2, -0.4, 0.6], Σ)
            x = rand(StableRNG(789990641), mvn)
            y = transport_to(StandardDist{Normal}(3), mvn)(x)
            @test transport_to(mvn, StandardDist{Normal}(3))(y) ≈ x
        end
    end
end
