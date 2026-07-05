# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsExt

using LinearAlgebra: Diagonal, diag, dot, cholesky

import Random
using Random: AbstractRNG, rand!

import DensityInterface
using DensityInterface: logdensityof, densityof

import MeasureBase
using MeasureBase: AbstractMeasure, AsMeasure, asmeasure
using MeasureBase: Lebesgue, Counting, ℝ
using MeasureBase: StdMeasure, StdUniform, StdExponential, StdLogistic, StdNormal
using MeasureBase: PowerMeasure, WeightedMeasure, SuperpositionMeasure, PushforwardMeasure
using MeasureBase: basemeasure, rootmeasure, testvalue, productmeasure, pushfwd, superpose
using MeasureBase: getdof, checked_arg, massof
using MeasureBase: transport_to, transport_def, transport_origin, from_origin, to_origin
using MeasureBase: NoTransportOrigin, NoTransport
using MeasureBase: Reshape
using MeasureBase: convert_realtype, firsttype, _fwddiff, @_adignore
import MeasureBase:
    _dist_params_numtype, _trafo_cdf_impl, _trafo_quantile_impl, _trafo_quantile_impl_generic
using MeasureBase: _pushfront, _pushback, _dropfront, _dropback, _rev_cumsum, _exp_cumsum_log

import Distributions
using Distributions: Distribution, VariateForm, ValueSupport, ContinuousDistribution
using Distributions: Univariate, Multivariate, ArrayLikeVariate, Continuous, Discrete
using Distributions: Uniform, Exponential, Logistic, Normal
using Distributions: MvNormal, AbstractMvNormal, Beta, Dirichlet
using Distributions: ReshapedDistribution, AbstractMixtureModel

import Statistics
import StatsBase
import StatsFuns
import PDMats

using IrrationalConstants: log2π, invsqrt2π

using HeterogeneousComputing: real_numtype

using Static: True, False, StaticInt, static, dynamic
using StaticThings: asnonstatic
using FillArrays: Fill, Ones, Zeros

using ArgCheck: @argcheck

using ArraysOfArrays: ArrayOfSimilarArrays, flatview

include("measure_interface.jl")
include("standard_dist.jl")
include("standard_uniform.jl")
include("standard_normal.jl")
include("distribution_measure.jl")
include("dist_vartransform.jl")
include("univariate.jl")
include("standardmv.jl")
include("product.jl")
include("reshaped.jl")
include("mixture.jl")
include("dirichlet.jl")
include("dirac.jl")

end # module MeasureBaseDistributionsExt
