
using MLStyle
using Random: AbstractRNG

export @measure

# A fold over ASTs. Example usage in `replace`
function foldast(leaf, branch; kwargs...)
    function go(ast)
        MLStyle.@match ast begin
            Expr(head, args...) => branch(head, map(go, args); kwargs...)
            x                   => leaf(x; kwargs...)
        end
    end

    return go
end

function replace(p, f, expr)
    leaf(s) = p(s) ? f(s) : s 
    branch(head, newargs) = Expr(head, newargs...)
    foldast(leaf, branch)(expr)
end

# from https://thautwarm.github.io/MLStyle.jl/latest/tutorials/capture.html
function capture(template, ex, action)
    let template = Expr(:quote, template)
        quote
            @match $ex begin 
                $template => $action
                _         => nothing
            end
        end 
    end
end

macro capture(template, ex, action)
    capture(template, ex, action) |> esc
end



function _measure(expr)
    @capture $μ($(p...)) expr begin
        q = quote
            struct $μ{N,T} <: ParameterizedMeasure{N}
                par :: NamedTuple{N,T}
            end

            $KeywordCalls.@kwstruct $expr
        end   
        
        if !isempty(p)
            # e.g. Normal(μ,σ) = Normal((μ=μ, σ=σ))
            pnames = QuoteNode.(p)
            push!(q.args, :($μ($(p...)) = $μ(NamedTuple{($(pnames...),)}(($(p...),)))))
        end
        
        return q
    end
end

"""
    @measure <declaration>
    
The <declaration> gives a measure and its default parameters, and specifies
its relation to its base measure. For example,

    @measure Normal(μ,σ)

declares the `Normal` is a measure with default parameters `μ and σ`. The result is equivalent to
```
struct Normal{N,T} <: ParameterizedMeasure{N}
    par :: NamedTuple{N,T}
end

KeywordCalls.@kwstruct Normal(μ,σ)

Normal(μ,σ) = Normal((μ=μ, σ=σ))
```

See [KeywordCalls.jl](https://github.com/cscherrer/KeywordCalls.jl) for details on `@kwstruct`.
"""
macro measure(expr)
    esc(_measure(expr))
end


# (@macroexpand @measure Normal(μ,σ) ≃ (1/sqrt2π) * Lebesgue(X)) |> MacroTools.prettify


using MLStyle

macro μσ_methods(ex)
    esc(_μσ_methods(ex))
end

function _μσ_methods(ex)
    @match ex begin
        :($dist($(args...))) => begin
            argnames = QuoteNode.(args)

            d_args = (:(d.$arg) for arg in args)
            quote

                function Base.rand(rng::AbstractRNG, T::Type, d::$dist{($(argnames...), :μ, :σ)})
                    d.σ * rand(rng, T, $dist($(d_args...))) + d.μ
                end

                function logdensity(d::$dist{($(argnames...), :μ, :σ)}, x)
                    z = (x - d.μ) / d.σ   
                    return logdensity($dist($(d_args...)), z) - log(d.σ)
                end

                function Base.rand(rng::AbstractRNG, T::Type, d::$dist{($(argnames...), :σ)})
                    d.σ * rand(rng, T, $dist($(d_args...)))
                end

                function logdensity(d::$dist{($(argnames...), :σ)}, x)
                    z = x / d.σ   
                    return logdensity($dist($(d_args...)), z) - log(d.σ) 
                end

                function Base.rand(rng::AbstractRNG, T::Type, d::$dist{($(argnames...), :μ)})
                    rand(rng, T, $dist($(d_args...))) + d.μ
                end

                function logdensity(d::$dist{($(argnames...), :μ)}, x)
                    z = x - d.μ
                    return logdensity($dist($(d_args...)), z)
                end
            end 
        end
    end
end

macro half(ex)
    esc(_half(ex))
end

function _half(ex)
    @match ex begin
        :($dist($(args...))) => begin
            halfdist = Symbol(:Half, dist)

            quote
                export $halfdist
                
                struct $halfdist{N,T} <: ParameterizedMeasure{N}
                    par :: NamedTuple{N,T}
                end
                
                unhalf(μ::$halfdist) = $dist(getfield(μ, :par))

                function basemeasure(μ::$halfdist) 
                    b = basemeasure(unhalf(μ))
                    lw = b.logweight
                    return WeightedMeasure(logtwo + lw, Lebesgue(ℝ₊))
                end
            
                function logdensity(μ::$halfdist, x)
                    return logdensity(unhalf(μ), x)
                end

                function Base.rand(rng::AbstractRNG, T::Type, μ::$halfdist)
                    return abs(rand(rng, T, unhalf(μ)))
                end

                (::$halfdist ≪ ::Lebesgue{ℝ₊}) = true
                (::Lebesgue{ℝ₊} ≪ ::$halfdist) = true
            end
        end
    end
end
