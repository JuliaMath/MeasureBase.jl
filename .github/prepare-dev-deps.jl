using Pkg

# Remove this source pin once HeterogeneousComputing 0.3 is registered.
withenv("JULIA_PKG_PRECOMPILE_AUTO" => "0") do
    Pkg.add(
        url = "https://github.com/oschulz/HeterogeneousComputing.jl",
        rev = "81c0392723ae2cdf4ae4a805dcc9d336752e9eff",
    )
end
