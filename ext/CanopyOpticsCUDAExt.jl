module CanopyOpticsCUDAExt

using CanopyOptics: FourSAILGeometrySet
using CUDA: CuArray

import CanopyOptics: _array_for_leaf, _geometry_for_leaf

_geometry_for_leaf(geoms::FourSAILGeometrySet{<:CuArray}, ::CuArray) = geoms

function _geometry_for_leaf(geoms::FourSAILGeometrySet, ::CuArray)
    return FourSAILGeometrySet(
        CuArray(geoms.μs),
        CuArray(geoms.μo),
        CuArray(geoms.dϕ),
        CuArray(geoms.dso),
        CuArray(geoms.ks),
        CuArray(geoms.ko),
        CuArray(geoms.sob),
        CuArray(geoms.sof),
        CuArray(geoms.sdb),
        CuArray(geoms.sdf),
        CuArray(geoms.dob),
        CuArray(geoms.dof),
        CuArray(geoms.ddb),
        CuArray(geoms.ddf),
    )
end

function _geometry_for_leaf(geoms::FourSAILGeometrySet{<:CuArray}, ::AbstractArray)
    return FourSAILGeometrySet(
        Array(geoms.μs),
        Array(geoms.μo),
        Array(geoms.dϕ),
        Array(geoms.dso),
        Array(geoms.ks),
        Array(geoms.ko),
        Array(geoms.sob),
        Array(geoms.sof),
        Array(geoms.sdb),
        Array(geoms.sdf),
        Array(geoms.dob),
        Array(geoms.dof),
        Array(geoms.ddb),
        Array(geoms.ddf),
    )
end

_array_for_leaf(x::CuArray, ::CuArray) = x
_array_for_leaf(x::AbstractArray, ::CuArray) = CuArray(x)
_array_for_leaf(x::CuArray, ::AbstractArray) = Array(x)

end # module CanopyOpticsCUDAExt
