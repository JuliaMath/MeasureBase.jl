# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

import Enzyme

@testset "Enzyme shape and probability derivatives" begin
    # The gamma shape derivative is an independent 50-digit reference
    # (also used in test_reference_values.jl); the beta derivatives follow
    # from I_x(a, 1) = x^a and I_x(1, b) = 1 - (1 - x)^b.
    a, x = 2.5, 2.75
    g = Gamma(a, 1)
    dα = -0.245741695389706
    gx = pdf(g, x)
    b = 0.3
    dβa, dβb = b * log(b), -(1 - b) * log1p(-b)
    gamma(a) = Gamma(a, 1; check_args = false)
    beta(a, b) = Beta(a, b; check_args = false)

    for (f, args, expected) in (
        ((a, x) -> smf(gamma(a), x), (a, x), (dα, gx)),
        ((a, x) -> csmf(gamma(a), x), (a, x), (-dα, -gx)),
        ((a, p) -> invsmf(gamma(a), p), (a, cdf(g, x)), (-dα / gx, 1 / gx)),
        ((a, q) -> invcsmf(gamma(a), q), (a, ccdf(g, x)), (-dα / gx, -1 / gx)),
        ((a, b, x) -> smf(beta(a, b), x), (1.0, 1.0, b), (dβa, dβb, 1.0)),
        ((a, b, p) -> invsmf(beta(a, b), p), (1.0, 1.0, b), (-dβa, -dβb, 1.0)),
    )
        reverse = Enzyme.gradient(Enzyme.Reverse, f, args...)
        @test all(isapprox.(reverse, expected; rtol = 1e-10))
        forward = only(
            Enzyme.autodiff(Enzyme.Forward, f, Enzyme.Duplicated.(args, one.(args))...),
        )
        @test forward ≈ sum(expected) rtol = 1e-10
    end

    # Constant shape arguments and zero adjoints at support boundaries.
    @test only(Enzyme.gradient(Enzyme.Reverse, p -> invsmf(g, p), cdf(g, x))) ≈ 1 / gx
    f(a) = smf(gamma(a), 0.0)
    @test only(Enzyme.gradient(Enzyme.Reverse, f, 2.5)) == 0
end

@testset "Enzyme beta constant shapes and Float32 tails" begin
    d, x = Beta(2.0, 3.0), 0.3
    for (f, z, expected) in
        ((z -> smf(d, z), x, pdf(d, x)), (p -> invsmf(d, p), cdf(d, x), inv(pdf(d, x))))
        reverse = only(Enzyme.gradient(Enzyme.Reverse, f, z))
        forward = only(Enzyme.autodiff(Enzyme.Forward, f, Enzyme.Duplicated(z, one(z))))
        @test reverse ≈ expected
        @test forward ≈ expected
    end

    # Q(p; a, 1) = p^(1/a); the mirrored upper-tail quantile uses
    # the complementary coordinate returned by the same inverse call.
    a, p = 50.0f0, 1.0f-6
    q = Float64(p)^(1 / Float64(a))
    dshape = -q * log(Float64(p)) / Float64(a)^2
    dp = q / (Float64(a) * Float64(p))
    for (f, expected) in (
        ((a, p) -> invsmf(Beta(a, 1.0f0; check_args = false), p), (dshape, dp)),
        ((a, p) -> invcsmf(Beta(1.0f0, a; check_args = false), p), (-dshape, -dp)),
    )
        args = (a, p)
        reverse = Enzyme.gradient(Enzyme.Reverse, f, args...)
        forward = only(
            Enzyme.autodiff(Enzyme.Forward, f, Enzyme.Duplicated.(args, one.(args))...),
        )
        @test all(isapprox.(reverse, expected; rtol = 64eps(Float32)))
        @test forward ≈ sum(expected) rtol = 64eps(Float32)
    end
end

@testset "Enzyme normal log-tail derivatives" begin
    for T in (Float32, Float64), x in T.((-16, 0, 16))
        for (f, reference) in ((logsmf, logcdf), (logcsmf, logccdf))
            expected =
                T(exp(logpdf(Normal(), Float64(x)) - reference(Normal(), Float64(x))))
            expected = f === logcsmf ? -expected : expected
            derivative = only(Enzyme.gradient(Enzyme.Reverse, z -> f(StdNormal(), z), x))
            @test derivative ≈ expected rtol = eps(T) * (4 + abs2(x))
        end
    end
end
