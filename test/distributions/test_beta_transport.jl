using Test, Distributions, MeasureBase
using MeasureBase: asmeasure, StdUniform, transport_to, smf, csmf, invsmf

@testset "Beta transport accuracy ($T)" for T in (Float32, Float64)
    p = T[1e-6, 0.01, 0.25, 0.5, 0.9]
    singular = asmeasure(Beta(T(1) / 2, T(1) / 2))
    expected = T.(sinpi.(BigFloat.(p) ./ 2) .^ 2)
    @test invsmf.(Ref(singular), p) ≈ expected rtol=64eps(T)

    polynomial = asmeasure(Beta(T(2), T(3)))
    threshold = T(3) / 7
    x = T[0.01, prevfloat(threshold), threshold, nextfloat(threshold), 0.9999]
    # Direct integration gives F=x²(6-8x+3x²), Q=(1-x)³(1+3x).
    xb = BigFloat.(x)
    @test smf.(Ref(polynomial), x) ≈ T.(xb .^ 2 .* (6 .- 8xb .+ 3xb .^ 2)) rtol=64eps(T)
    @test csmf.(Ref(polynomial), x) ≈ T.((1 .- xb) .^ 3 .* (1 .+ 3xb)) rtol=64eps(T)

    # A rounded complement must not erase an O(1) term multiplied by a large shape.
    b, x = T(1e7), T(1e-8)
    expected_small = T(-expm1(BigFloat(b) * log1p(-BigFloat(x)) + log1p(BigFloat(b) * x)))
    @test smf(asmeasure(Beta(T(2), b)), x) ≈ expected_small rtol=64eps(T)

    # The Beta coordinate rounds to one in Float32, but its complement is useful.
    ratio = asmeasure(Distributions.BetaPrime(T(1), T(0.02)))
    expected_ratio = T(expm1(-log1p(-BigFloat(T(0.3))) / BigFloat(T(0.02))))
    @test invsmf(ratio, T(0.3)) ≈ expected_ratio rtol=64eps(T)

    @test invsmf(asmeasure(Beta(T(1e-4), T(2))), T(0.3)) == 0
    @test invsmf(asmeasure(Beta(T(2), T(1e-4))), T(0.3)) == 1

    concentrated = asmeasure(Beta(T(1e5), T(1e5)))
    @test transport_to(concentrated, StdUniform())(T(0.5)) ≈ T(0.5) rtol=64eps(T)
end
