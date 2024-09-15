struct Reshape{SZA<:Tuple, SZB<:Tuple} <: Function
    new_shape::SZA
    old_shape::SZB
end

Base.@propagate_inbounds function (f::Reshape{SZA, SZB})(x::AbstractArray)
    @boundcheck @argcheck size(x) == f.old_shape
    return reshape(x, f.new_shape)
end

InverseFunctions.inverse(f::Reshape) = Reshape(f.old_shape, f.new_shape)

ChangesOfVariables.with_logabsdet_jacobian(f::Reshape, x::AbstractArray) = zero(real_numtype(typeof(x)))


