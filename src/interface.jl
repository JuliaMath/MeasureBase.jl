"""
    module MeasureBase.Interface

Test utilities for measure implementations. They are implemented in a
package extension, load `Test` (`using Test`) to use them.
"""
module Interface

using Reexport

@reexport using MeasureBase

using MeasureBase: basemeasure_depth, proxy
using MeasureBase: insupportof, basemeasure_sequence
using MeasureBase: transport_to

export test_interface
export test_transport
export test_smf
export basemeasure_depth
export proxy
export insupportof
export basemeasure_sequence

"""
    test_interface(μ)

Test that the base measure depth, test value, density and base measure of
`μ` are consistent and type stable, and that batched densities
(`logdensities`, densities of powers) agree with the point densities at
random variates of `μ` (at its test value if `μ` has infinite mass).
Requires `using Test`.
"""
function test_interface end

"""
    test_transport(ν, μ)

Test that `transport_to(ν, μ)` round-trips, is type stable and that its
log-abs-det Jacobian matches the densities of `μ` and `ν`, and that
batched transports (`transports`, broadcasting, `transports_with_ladj`)
agree with the point transport at random variates of `μ`. Requires
`using Test`.
"""
function test_transport end

"""
    test_smf(μ, n = 100)

Test the Stieltjes measure function `smf(μ, x)` and its inverse against
each other and against the measure of intervals. Requires `using Test`.
"""
function test_smf end

end # module Interface
