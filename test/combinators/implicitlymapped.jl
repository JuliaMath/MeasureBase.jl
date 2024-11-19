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

using Test
using MeasureBase
using Static: static
using Random: MersenneTwister

@testset "TakeAny" begin
    rng = MersenneTwister(42)
    
    @testset "Basic properties" begin
        take2 = TakeAny(2)
        take_static2 = TakeAny(static(2))
        
        # Test with various collection types
        arr = [1,2,3,4,5]
        @test length(take2(arr)) == 2
        @test length(take_static2(arr)) == 2
        
        # Test consistency
        @test take2(arr) == take2(arr)  # Same elements when called multiple times
        
        # Test with different sized inputs
        @test length(take2(1:10)) == 2
        @test length(take2(1:1)) == 1  # Should handle cases where input is smaller than n
    end

    @testset "Implicit mapping with TakeAny" begin
        # Create a kernel that produces a product measure
        kernel = par -> StdNormal()^3
        
        # Create mapped version that only looks at first two components
        mapped_kernel = ImplicitlyMapped(kernel, TakeAny(2))
        
        # Test with some parameter value
        par = 1.0
        full_measure = kernel(par)
        mapped_measure = explicit_kernel(mapped_kernel, rand(rng, full_measure))(par)
        
        # Check dimensions
        @test getdof(mapped_measure) == 2
        @test getdof(full_measure) == 3
        
        # Test consistency of mapping
        obs1 = rand(rng, full_measure)
        obs2 = rand(rng, full_measure)
        
        mapped1 = mapped_kernel.mapfunc(obs1)
        mapped2 = mapped_kernel.mapfunc(obs2)
        
        # Same elements should be selected consistently
        @test length(mapped1) == 2
        @test mapped_kernel.mapfunc(obs1) == mapped1  # Consistent mapping
    end
end