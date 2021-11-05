struct MapsTo{X,Y}
    x::X
    y::Y
end

↦(x::X, y::Y) where {X,Y} = MapsTo{X,Y}(x, y)

mapsto(x, y) = x ↦ y

Base.first(t::MapsTo) = t.x
Base.last(t::MapsTo) = t.y

Base.Pair(t::MapsTo) = t.x => t.y

Base.show(io::IO, t::MapsTo) = print(io, t.x, " ↦ ", t.y)

logdensity(d, t::MapsTo) = logdensity(d, t.y)
