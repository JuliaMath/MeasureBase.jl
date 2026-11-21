# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

# Checks of draws of discrete distributions against their probabilities.

import Distributions, Statistics

# A bound on the probability that draws of the discrete distribution `d`
# deviate from its probabilities as much as the draws `xs` do: the smallest
# two-sided binomial tail probability over the cells of their histogram,
# times the number of cells. The cells are up to 64 intervals of counts
# between the 1e-5 quantiles of `d` and the two tails beyond them.
function pmf_pvalue(xs::AbstractArray{<:Integer}, d::Distributions.Distribution)
    lo, hi = Int(Distributions.quantile(d, 1e-5)), Int(Distributions.quantile(d, 1 - 1e-5))
    w = cld(hi - lo + 1, 64)
    edges = lo:w:(hi + w)
    ncells = length(edges) - 1
    cellof(x) = x < lo ? 1 : x >= last(edges) ? ncells + 2 : (x - lo) ÷ w + 2
    counts = zeros(Int, ncells + 2)
    foreach(x -> counts[cellof(x)] += 1, xs)
    inner = [Distributions.cdf(d, edges[i+1] - 1) - Distributions.cdf(d, edges[i] - 1) for i in 1:ncells]
    probs = vcat(Distributions.cdf(d, lo - 1), inner, Distributions.ccdf(d, last(edges) - 1))
    tails = map(counts, probs) do c, p
        b = Distributions.Binomial(length(xs), clamp(p, 0, 1))
        2 * min(Distributions.cdf(b, c), Distributions.ccdf(b, c - 1))
    end
    return min(1, length(tails) * minimum(tails))
end

# The deviation of the mean of the draws `xs` from the mean of `d`, in
# standard deviations of the mean:
mean_zscore(xs::AbstractArray{<:Integer}, d::Distributions.Distribution) =
    (Statistics.mean(xs) - Statistics.mean(d)) / sqrt(Statistics.var(d) / length(xs))

# The deviation of the variance of the draws `xs` from the variance of `d`,
# in standard deviations of the sample variance:
var_zscore(xs::AbstractArray{<:Integer}, d::Distributions.Distribution) =
    (Statistics.var(xs) / Statistics.var(d) - 1) / sqrt((Distributions.kurtosis(d) + 2) / length(xs))
