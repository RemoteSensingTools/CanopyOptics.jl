"""
    _resolve_quadrature(quadrature::CanopyQuadrature, nQuad)

Return the explicit [`CanopyQuadrature`](@ref) unless the deprecated `nQuad`
keyword was supplied.

`nQuad` used to control both leaf-inclination and azimuth quadrature. It is
kept as a compatibility alias and maps to `CanopyQuadrature(nQuad, nQuad)`.
New code should pass `quadrature = CanopyQuadrature(...)`.
"""
function _resolve_quadrature(quadrature::CanopyQuadrature, nQuad)
    nQuad === nothing && return quadrature
    return CanopyQuadrature(nQuad, nQuad)
end
