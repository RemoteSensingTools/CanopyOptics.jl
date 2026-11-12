"""
    _resolve_quadrature(quadrature::CanopyQuadrature, nQuad)

Return `quadrature` when `nQuad === nothing`, or `CanopyQuadrature(nQuad, nQuad)`
when the compatibility keyword `nQuad` is supplied.
"""
function _resolve_quadrature(quadrature::CanopyQuadrature, nQuad)
    nQuad === nothing && return quadrature
    return CanopyQuadrature(nQuad, nQuad)
end
