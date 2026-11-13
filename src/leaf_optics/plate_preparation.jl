"""
    diffuse_plate_optics(material_depth, refractive_index)

Prepare one regular PROSPECT plate from finite material optical depth `k ≥ 0`
and refractive index `nr > 1`, using the same floating-point type. Return
`traversal_transmission`, `traversal_absorptance`, `effective_depth`,
`interface_entry_transmission` and `interface_exit_transmission`.

Traversal transmission is ``2 E_3(k)`` and absorption is evaluated directly near
zero. The effective exponential depth is its negative logarithm, computed from
direct loss near the conservative limit so a rounded unit transmission does not
erase weak absorption. This depth is an internal-source approximation, not
angle-resolved propagation. Both interfaces use hemispheric (``90^\\circ``) Fresnel
transmission, with ``t_{\\mathrm{exit}} = t_{\\mathrm{entry}} / n_r^2``. The entry convention differs from default
PROSPECT's first-surface 40-degree cone and is not a solar zenith angle.

This is host preparation, tested in Float32/Float64 over the native optical
database's refractive-index range. The nearly index-matched limit of the legacy
Fresnel expression is unsupported; invalid computed interface probabilities are
rejected. There is no CUDA or ForwardDiff contract for this helper. Use its
prepared values in [`RegularPlateEscape`](@ref) and [`ExponentialBirthProfile`](@ref)
for CPU/CUDA application. Zero traversal transmission from underflow is rejected
instead of floored; the scalar escape core can accept independently resolved
finite effective depths beyond this preparation domain.

No plate count, material composition, wavelength grid, intrinsic spectrum,
physiological yield or source normalization is chosen here. In particular, this
does not reconstruct a fractional PROSPECT remainder's internal interfaces.
"""
function diffuse_plate_optics(k::FT, nr::FT) where {FT<:AbstractFloat}
    isfinite(k) && k >= zero(FT) ||
        throw(ArgumentError("material optical depth must be finite and nonnegative"))
    isfinite(nr) && nr > one(FT) ||
        throw(ArgumentError("refractive index must be finite and greater than one"))
    attenuation = prospect_internal_attenuation(k)
    tau, loss = attenuation.transmission, attenuation.absorptance
    isfinite(tau) && zero(FT) < tau <= one(FT) || throw(
        DomainError(
            k,
            "unresolved traversal transmission; use higher precision or an independently resolved effective depth",
        ),
    )
    isfinite(loss) && zero(FT) <= loss <= one(FT) ||
        throw(DomainError(k, "invalid direct traversal absorption"))
    effective_depth = loss <= inv(FT(2)) ? -log1p(-loss) : -log(tau)
    entry = FT(calctav(FT(90), nr))
    exit = entry / (nr * nr)
    isfinite(entry) &&
    zero(FT) < entry <= one(FT) &&
    isfinite(exit) &&
    zero(FT) < exit <= one(FT) || throw(
        DomainError(
            nr,
            "invalid legacy Fresnel interface transmission; nearly matched indices are unsupported",
        ),
    )
    return (;
        traversal_transmission = tau,
        traversal_absorptance = loss,
        effective_depth,
        interface_entry_transmission = entry,
        interface_exit_transmission = exit,
    )
end
