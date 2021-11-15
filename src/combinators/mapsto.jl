struct MapsTo{X,Y}
    x::X
    y::Y
end

export ↦, mapsto

mapsto(x, y) = x ↦ y

↦(x::X, y::Y) where {X,Y} = MapsTo{X,Y}(x, y)

Base.first(t::MapsTo) = t.x
Base.last(t::MapsTo) = t.y

Base.Pair(t::MapsTo) = t.x => t.y

Base.show(io::IO, t::MapsTo) = print(t.x, " ↦ ", t.y)

logdensity(d, t::MapsTo) = logdensity(d, t.y)

Base.iterate(m::MapsTo) = (first(m), 1)
Base.iterate(m::MapsTo, ::Int) = (last(m), nothing)
Base.iterate(m::MapsTo, ::Nothing) = nothing