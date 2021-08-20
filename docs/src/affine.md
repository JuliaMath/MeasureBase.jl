# Affine Transformations

It's very common for measures to be parameterized by `μ` and `σ`, for example as in `Normal(μ=3,σ=4)` or `StudentT(ν=1,μ=3,σ=4)`. In this context, `μ` and `σ` do not always refer to the mean and standard deviation (the `StudentT` above is equivalent to a Cauchy, so both are undefined).

Rather, `μ` is a "location parameter", and `σ` is a "scale parameter". Together, these determine a transform 

```math
x → σx + μ
```

