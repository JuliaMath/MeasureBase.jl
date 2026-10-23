# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

module MeasureBaseDistributionsExt

using LinearAlgebra: Diagonal, Cholesky, LowerTriangular, UpperTriangular, Symmetric, diag, cholesky, logdet
import Adapt

using Random: AbstractRNG

import DensityInterface
using DensityInterface: logdensityof, densityof

import MeasureBase
using MeasureBase: AbstractMeasure, AsMeasure, asmeasure
using MeasureBase: Lebesgue, Counting, ℝ
using MeasureBase: StdMeasure, StdUniform, StdExponential, StdLogistic, StdNormal
using MeasureBase: _unit_interior
using MeasureBase: WeightedMeasure, SuperpositionMeasure, PushforwardMeasure
using MeasureBase: basemeasure, rootmeasure, testvalue, productmeasure, pushfwd, superpose
using MeasureBase: getdof, massof
using MeasureBase: transport_to, transport_def, transport_to_std, transport_from_std
using MeasureBase: batched_transport_to_std, batched_transport_from_std
using MeasureBase: Reshape
using MeasureBase: convert_realtype, _fwddiff, @_adignore
import MeasureBase:
    _dist_params_numtype, _trafo_logcdf_impl, _trafo_logccdf_impl,
    _trafo_quantile_impl, _trafo_cquantile_impl, _dist_quantile, _dist_cquantile
using MeasureBase: _pushfront, _pushback, _dropfront, _dropback, _rev_cumsum, _exp_cumsum_log
using MeasureBase: _gamma_cdf, _gamma_quantile, _beta_cdf, _beta_quantile, _gamma_logpdf, _beta_logpdf, _dualtag

import Distributions
using Distributions: Distribution, VariateForm, ValueSupport, ContinuousDistribution
using Distributions: Univariate, Multivariate, ArrayLikeVariate, Continuous, Discrete
using Distributions: Uniform, Exponential, Logistic, Normal
using Distributions: Cauchy, Laplace, LogNormal, Weibull, Gamma, Poisson, Bernoulli
using Distributions: MvNormal, AbstractMvNormal, Beta, Dirichlet
using Distributions: ReshapedDistribution, AbstractMixtureModel

import Statistics
import PDMats

using IrrationalConstants: log2π
using LogExpFunctions: logistic, log1pexp
using SpecialFunctions: loggamma, logbeta, gamma_inc, gamma_inc_inv, beta_inc, beta_inc_inv

using HeterogeneousComputing: real_numtype, GenContext, get_rng, get_precision

using Static: True, False, static
using StaticThings: SizeLike, asnonstatic
using FillArrays: Fill

using ArgCheck: @argcheck

using ArraysOfArrays: ArrayOfSimilarArrays, flatview

include("measure_interface.jl")
include("distribution_measure.jl")
include("dist_vartransform.jl")
include("univariate.jl")
include("families.jl")
include("standardmv.jl")
include("product.jl")
include("reshaped.jl")
include("mixture.jl")
include("dirichlet.jl")
include("multivariate.jl")
include("dirac.jl")
include("moments.jl")

end # module MeasureBaseDistributionsExt
