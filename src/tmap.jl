# A type-level map

function tmap(f::F, ::Type{T}) where {F,T<:Tuple}
    Tuple{map(f, T.types)...}
end

function tmap(f::F, ::Type{A}) where {F,T,A<:AbstractArray{T}}
    p = Tuple(A.parameters)
    C = constructorof(A)
    C{(f(T), Base.tail(p)...)...}
end

function tmap(f::F, ::Type{Base.Generator{G,I}}) where {F,G,I}
    Base.Generator{ComposedFunction{F,G}, I}
end

function tmap(f, ::Type{ProductMeasure{T}}) where {T}
    ProductMeasure{tmap(f, T)}
end
