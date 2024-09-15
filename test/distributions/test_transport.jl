# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

using DistributionMeasures
using Test

using LinearAlgebra
using InverseFunctions, ChangesOfVariables
using Distributions, ArraysOfArrays
import ForwardDiff, Zygote

using MeasureBase: transport_to, transport_def, transport_origin
using MeasureBase: StdExponential
using DistributionMeasures: _trafo_cdf, _trafo_quantile

include("getjacobian.jl")


@testset "test_distribution_transform" begin
    function test_back_and_forth(trg, src)
        @testset "transform $(typeof(trg).name) <-> $(typeof(src).name)" begin
            x = rand(src)
            y = transport_def(trg, src, x)
            src_v_reco = transport_def(src, trg, y)

            @test x ≈ src_v_reco
            
            f = x -> transport_def(trg, src, x)
            ref_ladj = logpdf(src, x) - logpdf(trg, y)
            @test ref_ladj ≈ logabsdet(getjacobian(f, x))[1]
        end
    end

    reshaped_rand(d::Distribution{Univariate}, n) = rand(d, n)
    reshaped_rand(d::Distribution{Multivariate}, n) = nestedview(rand(d, n))

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
        stduvnorm = StandardDist{Uniform}()

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

        @test isapprox(_trafo_cdf(Normal(Dual(0, 1, 0, 0), Dual(1, 0, 1, 0)), Dual(0.5, 0, 0, 1)), cdf(Normal(Dual(0, 1, 0, 0), Dual(1, 0, 1, 0)), Dual(0.5, 0, 0, 1)), rtol = 10^-6)
        @test isapprox(_trafo_cdf(Normal(0, 1), Dual(0.5, 1)), cdf(Normal(0, 1), Dual(0.5, 1)), rtol = 10^-6)

        @test isapprox(_trafo_quantile(Normal(0, 1), Dual(0.5, 1)), quantile(Normal(0, 1), Dual(0.5, 1)), rtol = 10^-6)
        @test isapprox(_trafo_quantile(Normal(Dual(0, 1, 0, 0), Dual(1, 0, 1, 0)), Dual(0.5, 0, 0, 1)), quantile(Normal(Dual(0, 1, 0, 0), Dual(1, 0, 1, 0)), Dual(0.5, 0, 0, 1)), rtol = 10^-6)
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
            (StandardNormal, StandardNormal)
            (Normal, StandardNormal)
            (StandardUniform, StandardUniform)
            (Uniform, StandardUniform)
        ]
            @test @inferred(transport_to(M, Weibull())) == transport_to(R(), Weibull())
            @test @inferred(transport_to(Weibull(), M)) == transport_to(Weibull(), R())
            @test @inferred(transport_to(M, MvNormal(float(I(5))))) == transport_to(R(5), MvNormal(float(I(5))))
            @test @inferred(transport_to(MvNormal(float(I(5))), M)) == transport_to(MvNormal(float(I(5))), R(5))
            @test @inferred(transport_to(M, StdExponential()^(2,3))) == transport_to(R(6), StdExponential()^(2,3))
            @test @inferred(transport_to(StdExponential()^(2,3), M)) == transport_to(StdExponential()^(2,3), R(6))
        end
    end
end
