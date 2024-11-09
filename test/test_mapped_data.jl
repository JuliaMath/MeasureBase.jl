using Test

using MeasureBase

using StaticArrays: SVector
using Static: static
using AffineMaps, PropertyFunctions

@testset "implicitlymapped" begin
    @testset "TakeAny" begin
        V = [3, 2, 4, 2, 7, 5, 6]
        mV = [3, 2, 4, 2]
        S = Set(V)
        mS = Set(mV)
        SV = SVector(V...)
        mSV = SVector(mV...)

        @test @inferred(MeasureBase.TakeAny(4)(V)) == mV
        @test @inferred(MeasureBase.TakeAny(static(4))(V)) == mV
        tS = @inferred(MeasureBase.TakeAny(4)(S))
        @test tS isa Set && length(tS) == 4 && all(x -> x in S, tS)
        @test @inferred(MeasureBase.TakeAny(static(4))(S)) == MeasureBase.TakeAny(4)(S)
        @test @inferred(MeasureBase.TakeAny(static(4))(SV)) === mSV
        @test @inferred(MeasureBase.TakeAny(4)(SV)) == mV
        @test @inferred(MeasureBase.TakeAny(4)(V)) == mV
    end

    function test_implicitly_mapped(
        label,
        f_kernel,
        ref_mapfunc,
        ref_mappedkernel,
        par,
        orig_obs,
        obs,
    )
        @testset "$label" begin
            im_measure = @inferred Marginalized(f_kernel(par))
            im_kernel = @inferred Marginalized(f_kernel)
            mapfunc = @inferred explicit_mapfunc(im_measure, obs)
            mapped_measure = @inferred explicit_measure(im_measure, obs)
            mapped_likelihood = @inferred Likelihood(im_kernel, obs)

            @test mapfunc == ref_mapfunc
            @test @inferred(mapfunc(orig_obs)) == obs
            @test mapped_measure == ref_mappedkernel(par)

            @test @inferred(logdensityof(im_measure, obs)) ≈
                  logdensityof(mapped_measure, obs)
            @test @inferred(logdensityof(mapped_likelihood, par)) ≈
                  logdensityof(Likelihood(ref_mappedkernel, obs), par)
        end
    end

    f_kernel =
        par -> productmeasure(
            map(
                m -> pushfwd(Mul(par), m),
                (a = StdUniform(), b = StdNormal(), c = StdExponential()),
            ),
        )
    ref_mapfunc = @pf (; $a, $c)
    ref_mappedkernel =
        par -> productmeasure(
            map(m -> pushfwd(Mul(par), m), (a = StdUniform(), c = StdExponential())),
        )
    par = 4.2
    orig_obs = (a = 0.7, b = 2.1, c = 1.2)
    obs = (a = 0.7, c = 1.2)
    test_implicitly_mapped(
        "marginalized nt",
        f_kernel,
        ref_mapfunc,
        ref_mappedkernel,
        par,
        orig_obs,
        obs,
    )

    f_kernel = par -> pushfwd(Mul(par), StdNormal())^7
    ref_mapfunc = MeasureBase.TakeAny(3)
    ref_mappedkernel = par -> pushfwd(Mul(par), StdNormal())^3
    par = 4.2
    orig_obs = [9.4, -7.3, 1.0, -2.9, 1.9, 4.7, 0.5]
    obs = [9.4, -7.3, 1.0]
    test_implicitly_mapped(
        "marginalized nt",
        f_kernel,
        ref_mapfunc,
        ref_mappedkernel,
        par,
        orig_obs,
        obs,
    )
end
