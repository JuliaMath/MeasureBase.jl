_LATENT_DOCSTRING = """
Some Probabilistic Programming Languages (PPLs) like Tilde.jl make a distinction
between a _latent space_, often a namespace represented as a named tuple, and
the space containing the return value, which we refer to as the _maifest space_.
The distinction is that computations are done in terms of the latent space,
while the resulting value is in the manifest space.

To simplify many manipulations involving these concepts, we introduce the
concept of a _joint space_. For example, suppose `m()` is a measure with latent
space `NamedTuple{(:a, :b)}` that returns `a - b`, so the latent value `(a = 3,
b = 4)` is mapped to the manifest value `0.75`. Then the corresponding value in
the joint space is the pair `(a = 3, b = 4) => 0.75`.

One of the many goals of probabilistic programming is to blur the line between
"built in" measures like `Normal()` and those defined in terms of a model from a
PPL. To accommodate this, we extend these concepts to general measures. 

For many measures, it's convenient to work directly in the manifest space, and
there's no need for such separation. However, it's important to be able to
manipulate measures programmatically, with minimal special cases. Because of
this, we introduce fall-back methods

    latentof(m) = m
    manifestof(m) = m

The default implementation of `jointof` is then a push-forward through the
function `x -> (x => x)`. For example,

    julia> rand(MeasureBase.jointof(StdUniform()))
    0.346439=>0.346439    
"""

"""
    latentof(m)

$_LATENT_DOCSTRING
"""
latentof(m) = m

"""
    manifestof(m)

$_LATENT_DOCSTRING
"""
manifestof(m) = m

"""
    jointof(m)

$_LATENT_DOCSTRING
"""
function jointof(m)
    fwd(x) = x => x
    
    function back(p::Pair)
        x,y = p
        @assert x === y
        return x
    end

    PushforwardMeasure(fwd, back, m, NoVolCorr())
end
