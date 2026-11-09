using Test

using MeasureBase
using MeasureBase: pushfwd, StdUniform, StdExponential, StdLogistic
using MeasureBase: pushfwd, PushforwardMeasure
using MeasureBase: transport_to, unsafe_logdensityof
using MeasureBase: productmeasure, mbind, superpose
import AffineMaps
import Zygote
using InverseFunctions: inverse
using ChangesOfVariables: with_logabsdet_jacobian
using Statistics: var
using DensityInterface: logdensityof
using LogExpFunctions
using SpecialFunctions: erfc, erfcinv
import InverseFunctions: inverse, setinverse
using IrrationalConstants: invsqrt2, sqrt2
using IntervalSets
import ChangesOfVariables: with_logabsdet_jacobian
using MeasureBase.Interface: transport_to, test_transport
using ArraysOfArrays: flatview
using Random: Xoshiro

Φ(z) = erfc(-z * invsqrt2) / 2
Φinv(p) = -erfcinv(2 * p) * sqrt2

with_logabsdet_jacobian(::typeof(Φ), z) = (Φ(z), logdensityof(StdNormal(), z))

function with_logabsdet_jacobian(::typeof(Φinv), p)
    z = Φinv(p)
    (z, -logdensityof(StdNormal(), z))
end

inverse(::typeof(Φ)) = Φinv
inverse(::typeof(Φinv)) = Φ

function test_pushfwd(f, μ, ν_ref)
    @testset "pushfwd($f, $μ)" begin
        @inferred(pushfwd(f, μ))
        ν = pushfwd(f, μ)

        test_transport(ν, ν_ref)
        test_transport(ν_ref, ν)

        y = rand(ν_ref)
        @test isapprox(@inferred(logdensityof(ν, y)), logdensityof(ν_ref, y), atol = 1e-10)
        @test isapprox(
            @inferred(unsafe_logdensityof(ν, y)),
            unsafe_logdensityof(ν_ref, y),
            atol = 1e-10,
        )

        @test isapprox(var(rand(ν^(10^5))), var(ν_ref), rtol = 0.05)

        @test abs(transport_to(StdLogistic(), ν)(y)) ≈
              abs(transport_to(StdLogistic(), ν_ref)(y))
    end
end

@testset "transformedmeasure.jl" begin
    # (f, μ, ν_ref), so that pushfwd(f, μ) ≅ ν_ref
    triples = [
        ((-) ∘ log, StdUniform(), StdExponential())
        (exp ∘ (-), StdExponential(), StdUniform())
        (logit, StdUniform(), StdLogistic())
        (logistic, StdLogistic(), StdUniform())
        (Φ, StdNormal(), StdUniform())
        (Φinv, StdUniform(), StdNormal())
    ]

    for (f, μ, ν_ref) in triples
        test_pushfwd(f, μ, ν_ref)
        test_pushfwd(setinverse(f, inverse(f)), μ, ν_ref)
    end

    @testset "Pushforward-of-pushforward" begin
        for (f, μ, ν_ref) in triples
            finv = inverse(f)
            test_pushfwd(finv, pushfwd(f, μ), μ)
        end
    end
end

using IrrationalConstants: loghalf

@testset "Half" begin
    μ = Half(StdNormal())

    # Test logdensityof for positive values
    @test logdensityof(μ, 1.0) ≈ logdensityof(StdNormal(), 1.0) - loghalf

    # Test logdensityof for negative values (should return -Inf)
    @test logdensityof(μ, -1.0) == -Inf

    # Test logdensityof for zero
    @test logdensityof(μ, 0.0) ≈ logdensityof(StdNormal(), 0.0) - loghalf

    @test logdensities(μ, [-1.0, 1.0]) ≈ [-Inf, logdensityof(μ, 1.0)]
    @test MeasureBase.insupportof(μ, -1.0) == false && MeasureBase.insupportof(μ, 1.0) == true
    @test logdensityof(Half(2 * StdNormal()), 1.0) ≈ log(2) + logdensityof(μ, 1.0)

    # Half requires scalar variates and parents symmetric about zero:
    @test_throws ArgumentError Half(StdExponential())
    @test_throws ArgumentError Half(StdUniform())
    @test_throws ArgumentError Half(StdNormal()^3)

    # The fold would double the mass of atoms at zero, parents with atoms
    # are refused, continuous mixtures fold:
    @test_throws ArgumentError Half(Dirac(0.0))
    @test_throws ArgumentError Half(superpose(StdNormal(), Dirac(0.0)))
    φ(x) = exp(logdensityof(StdNormal(), x))
    mix = superpose(0.5 * pushfwd(AffineMaps.Add(-1.0), StdNormal()), 0.5 * pushfwd(AffineMaps.Add(1.0), StdNormal()))
    @test logdensityof(Half(mix), 0.7) ≈ log(φ(1.7) + φ(-0.3))

    # Half of a weighted parent has its mass, its transports map the
    # normalized law (Φ⁻¹(3/4) and 2Φ(1/2) - 1), also within truncations:
    h = Half(2 * StdNormal())
    @test massof(h) ≈ 2 && smf(h, 0.5) ≈ 2 * (2 * Φ(0.5) - 1) && smf(h, -1.0) == 0
    @test transport_to(h, StdUniform())(0.5) ≈ Φinv(0.75)
    @test transport_to(StdUniform(), h)(0.5) ≈ 2 * Φ(0.5) - 1
    @test transport_to(h, StdUniform())(0.5f0) isa Float32
    test_transport(h, StdUniform())
    for (a, b) in ((0.2, 1.5), (3.0, Inf))
        t, t1 = truncate(h, a .. b), truncate(Half(StdNormal()), a .. b)
        @test massof(t) ≈ 2 * massof(t1)
        @test transport_to(t, StdUniform())(0.3) ≈ transport_to(t1, StdUniform())(0.3)
        @test transport_to(StdUniform(), t)(a + 0.1) ≈ transport_to(StdUniform(), t1)(a + 0.1)
    end
end

using Test
using MeasureBase
using MeasureBase: gettransform
using InverseFunctions
using StaticArrays
using ChangesOfVariables
using ChangesOfVariables: setladj

@testset "PushforwardMeasure" begin
    # Test basic pushforward construction
    μ = StdNormal()
    f = exp
    ν = pushfwd(f, μ)

    @test ν isa PushforwardMeasure
    @test MeasureBase.paramsof(ν) == NamedTuple() && MeasureBase.paramnames(ν) == ()

    # TODO: How do we access the parent?
    # @test parent(ν) === μ
    @test gettransform(ν) === f

    # Test logdensity with AdaptRootMeasure
    x = 0.5
    y = f(x)
    ld = logdensityof(ν, y)
    @test !isnan(ld)
    @test isfinite(ld)

    # Test logdensity with PushfwdRootMeasure
    ν_no_corr = pushfwd(f, μ, PushfwdRootMeasure())
    ld_no_corr = logdensityof(ν_no_corr, y)
    @test !isnan(ld_no_corr)
    @test isfinite(ld_no_corr)

    # Test non-bijective pushforward
    g = x -> x^2
    @test MeasureBase.getdof(pushfwd(g, μ)) isa MeasureBase.NoDOF
    @test_throws ArgumentError logdensityof(pushfwd(g, μ), 1.0)
    @test_throws ArgumentError logdensityof(pushfwd(g, μ, PushfwdRootMeasure()), 1.0)

    # Test edge cases for _combine_logd_with_ladj
    @test MeasureBase._combine_logd_with_ladj(-Inf, Inf) == -Inf  # Zero density wins
    @test MeasureBase._combine_logd_with_ladj(1.0, -Inf) == -Inf
    @test logdensityof(pushfwd(log, StdExponential()), -Inf) == -Inf

    # Test composition of pushforwards
    h = x -> exp(x)
    ν_comp = pushfwd(h, ν)
    @test ν_comp isa PushforwardMeasure

    # TODO: How do we access the parent?
    # @test parent(ν_comp) === μ

    # Test identity pushforward
    ν_id = pushfwd(identity, μ)
    @test ν_id === μ

    # Test rootmeasure behavior
    @test rootmeasure(ν) === rootmeasure(μ)  # AdaptRootMeasure
    @test rootmeasure(ν_no_corr) isa PushforwardMeasure  # PushfwdRootMeasure

    # Pushforwards adapting the root have the origin's root as their base,
    # pushforwards of the root pull weights out of the pushed-forward base:
    @test basemeasure(ν) === rootmeasure(μ)
    @test basemeasure(ν_no_corr) isa WeightedMeasure{<:Any,<:PushforwardMeasure}
    @test basemeasure(ν_no_corr).base.style isa PushfwdRootMeasure

    # Both styles describe the same measure:
    @test logdensity_rel(ν, ν_no_corr, 2.0) ≈ 0 atol = 1e-12
    @test logdensity_rel(ν_no_corr, ν, 2.0) ≈ 0 atol = 1e-12
    @test logdensity_rel(ν, rootmeasure(μ), 2.0) ≈ logdensityof(μ, log(2.0)) - log(2.0)
    E = Base.BroadcastFunction(exp)
    νv, νv_no_corr = pushfwd(E, μ^2), pushfwd(E, μ^2, PushfwdRootMeasure())
    @test logdensity_rel(νv, νv_no_corr, [2.0, 3.0]) ≈ 0 atol = 1e-12
    @test logdensities(νv, [2.0 1.0; 3.0 0.5]) ≈ [logdensityof(νv, [2.0, 3.0]), logdensityof(νv, [1.0, 0.5])]

    # Rank-changing pushforwards reshape the root:
    mr = MeasureBase.mreshape(μ^6, (2, 3))
    @test basemeasure(mr) == MeasureBase.LebesgueBase()^(2, 3)
    xr = reshape(collect(0.1:0.1:0.6), 2, 3)
    @test logdensityof(mr, xr) ≈ logdensityof(μ^6, vec(xr))

    # Mixtures of pushforward location families:
    shifted = pushfwd(setinverse(setladj(x -> x + 1, zero), setladj(y -> y - 1, zero)), μ)
    s = superpose([μ, shifted])
    @test logdensityof(s, 0.3) ≈ log(exp(logdensityof(μ, 0.3)) + exp(logdensityof(μ, -0.7)))
    @test sum(x -> exp(logdensityof(s, x)), -12:0.01:13) * 0.01 ≈ 2

    # Test massof
    # TODO: mass interface is very incomplete
    # @test massof(ν) == massof(μ)

    # Test rand
    @test rand(ν) isa Real
    @test insupportof(ν, rand(ν)) != false

    # Test pullback
    pb = pullbck(f, ν)
    @test pb isa PushforwardMeasure
    @test logdensityof(pb, y) ≈ logdensityof(μ, y)

    # Test deprecated pullback
    @test_deprecated MeasureBase.pullback(f, μ)
end

@testset "pushforwards of counting measures" begin
    # Injective maps keep the masses of atoms, without a volume element:
    D = MeasureBase.Dirac
    a = superpose(2D(1.0), 3D(2.0))
    p = pushfwd(AffineMaps.Mul(2.0), a)
    @test p == superpose(2D(2.0), 3D(4.0))
    @test logdensityof(p, 2.0) ≈ log(2) && logdensityof(p, 4.0) ≈ log(3)
    @test logdensityof(p, 3.0) == -Inf
    @test sum(y -> exp(logdensityof(p, y)), 0.0:0.5:6.0) ≈ massof(a)
    @test pushfwd(AffineMaps.Mul(2.0), Counting(ℤ)).style isa PushfwdRootMeasure
    @test pushfwd(AffineMaps.Mul(2.0), StdNormal()).style isa AdaptRootMeasure
end

@testset "PushFwdStyle types" begin
    @test AdaptRootMeasure() isa PushFwdStyle
    @test PushfwdRootMeasure() isa PushFwdStyle
    @test MeasureBase.WithVolCorr === AdaptRootMeasure
    @test MeasureBase.NoVolCorr === PushfwdRootMeasure

    @testset "output size of pushforwards of tuple products" begin
        Pt = productmeasure((StdNormal(), StdUniform()^2))
        ν = pushfwd(x -> vcat(x[1], x[2]), Pt)
        @test MeasureBase.mspace_flatsize(ν) == (3,)
        @test MeasureBase.mspace_ndims(ν) == 1
        @test MeasureBase.mspace_flatsize(pushfwd(x -> x, mbind(x -> StdNormal()^(x > 0 ? 1 : 2), StdNormal()))) isa MeasureBase.NoMSpaceElementSize
    end

    @testset "construction inside differentiated functions" begin
        flat(x::Tuple) = vcat(x[1], x[2])
        unflat(v::AbstractVector) = (v[1], v[2:end])
        InverseFunctions.inverse(::typeof(flat)) = unflat
        InverseFunctions.inverse(::typeof(unflat)) = flat
        ChangesOfVariables.with_logabsdet_jacobian(::typeof(flat), x) = (flat(x), zero(eltype(x[2])))
        ChangesOfVariables.with_logabsdet_jacobian(::typeof(unflat), v) = (unflat(v), zero(eltype(v)))
        Pt = productmeasure((StdNormal(), StdUniform()^2))
        g(v) = logdensityof(pushfwd(flat, Pt), v)
        v = [0.3, 0.2, 0.7]
        @test Zygote.gradient(g, v)[1] ≈ [-0.3, 0.0, 0.0]
    end
end


# Doubling that counts its evaluations:
const _doubling_calls = Ref(0)
struct CountedDoubling end
struct CountedHalving end
(::CountedDoubling)(x) = (_doubling_calls[] += 1; 2 .* x)
(::CountedHalving)(y) = y ./ 2
InverseFunctions.inverse(::CountedDoubling) = CountedHalving()
InverseFunctions.inverse(::CountedHalving) = CountedDoubling()
ChangesOfVariables.with_logabsdet_jacobian(f::CountedHalving, y) = (f(y), -length(y) * log(2))

@testset "variate sizes of pushforwards" begin
    using MeasureBase: mspace_ndims, mspace_flatsize, fixed_stream_size, restrict, logdensities, StdNormal
    using AffineMaps: Mul
    using StaticArrays: SMatrix, Size
    using Static: static

    # Origins without random variates: construction and point densities
    # need no variate size, batches take it from the type `f` returns or
    # from a declaration:
    r = restrict(x -> x > 0, StdNormal())
    ν = pushfwd(Mul(2.0), r)
    @test logdensityof(ν, 1.0) ≈ logdensityof(StdNormal(), 0.5) - log(2)
    @test logdensityof(ν, -1.0) == -Inf
    @test @inferred(mspace_ndims(typeof(ν))) == 0
    @test logdensities(ν, [1.0, -1.0, 3.0]) ≈ [logdensityof(ν, y) for y in (1.0, -1.0, 3.0)]
    r2 = restrict(x -> all(>(0), x), StdNormal()^2)
    ν2 = pushfwd(CountedDoubling(), r2)
    @test logdensityof(ν2, [1.0, 2.0]) ≈ logdensityof(StdNormal()^2, [0.5, 1.0]) - 2 * log(2)
    @test_throws ArgumentError logdensities(ν2, [1.0 2.0; 2.0 1.0])
    ν2d = pushfwd(CountedDoubling(), r2; varsize = (2,))
    @test @inferred(mspace_ndims(typeof(ν2d))) == 1 && fixed_stream_size(typeof(ν2d)) == true
    @test logdensities(ν2d, [1.0 2.0; 2.0 -1.0]) ≈ [logdensityof(ν2, [1.0, 2.0]), -Inf]
    @test_throws ArgumentError pushfwd(Mul(2.0), r; varsize = 3)

    # Sizes are learned on demand, the function is not evaluated at
    # construction or for point densities:
    _doubling_calls[] = 0
    μ3 = StdNormal()^3
    ν3 = pushfwd(CountedDoubling(), μ3)
    y = [0.2, 0.4, 0.6]
    @test logdensityof(ν3, y) ≈ logdensityof(μ3, y ./ 2) - 3 * log(2)
    @test _doubling_calls[] == 0
    @test @inferred(mspace_ndims(typeof(ν3))) == 1
    @test @inferred(mspace_flatsize(ν3)) == (3,)
    @test _doubling_calls[] == 1
    Y = [0.2 0.4; 0.4 -0.8; 0.6 1.0]
    @test logdensities(ν3, Y) ≈ [logdensityof(ν3, Y[:, j]) for j in 1:2]
    @test_throws ArgumentError logdensities(ν3, Y[1:2, :])

    # Static sizes follow from the type `f` returns for statically sized
    # variates:
    n_calls = _doubling_calls[]
    A = SMatrix{2,2}(1.0, 0.5, 0.0, 2.0)
    νs = pushfwd(Mul(A), StdNormal()^static(2))
    @test @inferred(mspace_flatsize(typeof(νs))) === Size(2)
    @test @inferred(mspace_flatsize(typeof(pushfwd(CountedDoubling(), μ3; varsize = (static(3),))))) === Size(3)
    @test @inferred(mspace_flatsize(typeof(pushfwd(CountedDoubling(), StdNormal()^static(3))))) === Size(3)
    @test _doubling_calls[] == n_calls

    # Origins with nested variates: sizes are learned from the image of the
    # test value, functions that refuse lazy arrays need a declared size:
    o = (StdNormal()^2)^3
    @test mspace_flatsize(pushfwd(x -> reduce(vcat, x), o)) == (6,)
    dense_only(x) = x isa Vector ? reduce(vcat, x) : throw(ArgumentError("dense variates only"))
    @test_throws "varsize" mspace_flatsize(pushfwd(dense_only, o))
    @test mspace_flatsize(pushfwd(dense_only, o; varsize = (6,))) == (6,)
end

# A rotation of the plane with a fused batched implementation, counting
# the calls of the single-variate form:
const _rotation_calls = Ref(0)
struct Rotation2D
    θ::Float64
end
_rotmatrix(f::Rotation2D) = [cos(f.θ) -sin(f.θ); sin(f.θ) cos(f.θ)]
(f::Rotation2D)(x) = (_rotation_calls[] += 1; _rotmatrix(f) * x)
InverseFunctions.inverse(f::Rotation2D) = Rotation2D(-f.θ)
ChangesOfVariables.with_logabsdet_jacobian(f::Rotation2D, x) = (f(x), 0.0)
MeasureBase.batched_apply(f::Rotation2D, X) = _rotmatrix(f) * flatview(X)

@testset "variates passed to pushforward functions" begin
    # Functions receive the variates as `rand` returns them, single variates
    # and batches alike:
    o = (StdNormal()^2)^3
    s = pushfwd(x -> sum(x), o; varsize = (2,))
    x = rand(o)
    @test MeasureBase.batched_apply(x -> sum(x), [x, x]) ≈ [sum(x) sum(x)]
    @test rand(Xoshiro(1), s) ≈ sum(rand(Xoshiro(1), o))
    Y = rand(Xoshiro(2), s^5)
    @test size(Y) == (5,) && Y ≈ [sum(x) for x in rand(Xoshiro(2), o^5)]

    # Elementwise functions apply to the numbers of nested variates, which
    # keep their layout:
    p = pushfwd(Base.BroadcastFunction(exp), o)
    y = rand(Xoshiro(3), p)
    @test y ≈ [exp.(v) for v in rand(Xoshiro(3), o)]
    @test rand(Xoshiro(4), p^4) ≈ [[exp.(v) for v in x] for x in rand(Xoshiro(4), o^4)]
    @test logdensityof(p, y) ≈ logdensityof(o, [log.(v) for v in y]) - sum(v -> sum(log, v), y)
    Yp = rand(p^4)
    @test size(Yp) == (4,) && all(v -> length(v) == 3 && length(v[1]) == 2, Yp)
    @test logdensities(p, Yp) ≈ [logdensityof(p, Yp[i]) for i in 1:4]
    f = transport_to(StdUniform()^6, p)
    @test inverse(f)(f(y)) ≈ y
    U = rand((StdUniform()^6)^3)
    @test MeasureBase.transports(f, MeasureBase.transports(inverse(f), U)) ≈ U
    @test MeasureBase.batched_apply(Base.BroadcastFunction(exp), [[0.0, 1.0], [2.0, 3.0]]) ≈ exp.([0.0 2.0; 1.0 3.0])

    # Specializations of `batched_apply` are used for batches:
    r = Rotation2D(0.3)
    ν = pushfwd(r, StdNormal()^2; varsize = (2,))
    _rotation_calls[] = 0
    V = rand(ν^1000)
    @test _rotation_calls[] == 0
    @test size(flatview(V)) == (2, 1000)
    @test isapprox(sum(flatview(V); dims = 2) / 1000, zeros(2), atol = 0.15)
    Z = MeasureBase.transports(transport_to(StdNormal()^2, ν), V)
    @test Z ≈ [_rotmatrix(inverse(r)) * v for v in V]
    @test _rotation_calls[] == 0
    @test rand(ν) isa Vector{Float64} && _rotation_calls[] == 1
end
