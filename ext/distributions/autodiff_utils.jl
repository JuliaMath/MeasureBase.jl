# This file is a part of MeasureBase.jl, licensed under the MIT License (MIT).

@inline _adignore_call(f) = f()
@inline _adignore_call_pullback(@nospecialize ΔΩ) = (NoTangent(), NoTangent())
ChainRulesCore.rrule(::typeof(_adignore_call), f) = _adignore_call(f), _adignore_call_pullback

macro _adignore(expr)
    :(_adignore_call(() -> $(esc(expr))))
end


function _pushfront(v::AbstractVector, x)
    T = promote_type(eltype(v), typeof(x))
    r = similar(v, T, length(eachindex(v)) + 1)
    r[firstindex(r)] = x
    r[firstindex(r)+1:lastindex(r)] = v
    r
end

function ChainRulesCore.rrule(::typeof(_pushfront), v::AbstractVector, x)
    result = _pushfront(v, x)
    function _pushfront_pullback(thunked_ΔΩ)
        ΔΩ = unthunk(thunked_ΔΩ)
        (NoTangent(), ΔΩ[firstindex(ΔΩ)+1:lastindex(ΔΩ)], ΔΩ[firstindex(ΔΩ)])
    end
    return result, _pushfront_pullback
end


function _pushback(v::AbstractVector, x)
    T = promote_type(eltype(v), typeof(x))
    r = similar(v, T, length(eachindex(v)) + 1)
    r[lastindex(r)] = x
    r[firstindex(r):lastindex(r)-1] = v
    r
end

function ChainRulesCore.rrule(::typeof(_pushback), v::AbstractVector, x)
    result = _pushback(v, x)
    function _pushback_pullback(thunked_ΔΩ)
        ΔΩ = unthunk(thunked_ΔΩ)
        (NoTangent(), ΔΩ[firstindex(ΔΩ):lastindex(ΔΩ)-1], ΔΩ[lastindex(ΔΩ)])
    end
    return result, _pushback_pullback
end


_dropfront(v::AbstractVector) = v[firstindex(v)+1:lastindex(v)]

_dropback(v::AbstractVector) = v[firstindex(v):lastindex(v)-1]


_rev_cumsum(xs::AbstractVector) = reverse(cumsum(reverse(xs)))

function ChainRulesCore.rrule(::typeof(_rev_cumsum), xs::AbstractVector)
    result = _rev_cumsum(xs)
    function _rev_cumsum_pullback(ΔΩ)
        ∂xs = @thunk cumsum(unthunk(ΔΩ))
        (NoTangent(), ∂xs)
    end
    return result, _rev_cumsum_pullback
end


# Equivalent to `cumprod(xs)``:
_exp_cumsum_log(xs::AbstractVector) = exp.(cumsum(log.(xs)))

function ChainRulesCore.rrule(::typeof(_exp_cumsum_log), xs::AbstractVector)
    result = _exp_cumsum_log(xs)
    function _exp_cumsum_log_pullback(ΔΩ)
        ∂xs = inv.(xs) .* _rev_cumsum(exp.(cumsum(log.(xs))) .* unthunk(ΔΩ))
        (NoTangent(), ∂xs)
    end
    return result, _exp_cumsum_log_pullback
end
