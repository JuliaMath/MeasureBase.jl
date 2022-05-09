struct Splat{F}
    f::F
end

function (s::Splat{F})(x) where {F}
    s.f(x...)
end

unsplat(s::Splat) = s.f

splat(f) = Splat(f)
