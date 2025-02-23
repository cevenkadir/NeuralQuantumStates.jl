using StaticArrays: MArray

function get_array_type_params(
    array_type::Type{MArray}, T::Type, S::NTuple{D,Integer}
) where {D}
    return Tuple{S...}, T, D, prod(S)
end



function get_array_type_params(
    array_type::Type{Array}, T::Type, S::NTuple{D,Integer}
) where {D}
    return T, D
end
