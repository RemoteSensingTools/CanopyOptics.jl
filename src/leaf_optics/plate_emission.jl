"""
    AbstractLeafBirthProfile

Normalized internal photon birth-depth distribution, independent of yield and spectrum.
"""
abstract type AbstractLeafBirthProfile end

"""
    UniformBirthProfile()

Uniform photon production over dimensionless depth `0 ≤ z ≤ 1`.
"""
struct UniformBirthProfile <: AbstractLeafBirthProfile end

@doc raw"""
    ExponentialBirthProfile(excitation_depth, forward_amplitude, backward_amplitude)

Conditional birth density proportional to
``f\,\exp(-a z) + g\,\exp\bigl(-a(1 - z)\bigr)``, where ``a`` is `excitation_depth`,
``(f, g)`` are `forward_amplitude` and `backward_amplitude`, and ``z = 0`` is the
fixed anatomical front face.
These amplitudes weight excitation-depth contributions, not emission directions.
The host constructor conditions finite nonnegative relative amplitudes on their
positive sum, using scaling to avoid overflow and retaining both proportions
when representable (extreme amplitude ratios can underflow).
Depth must be finite and nonnegative; all three arguments share the same `FT`.

For a homogeneous, front-illuminated regular plate, the relative amplitudes are
``(1,\, r_{\mathrm{exc}}\exp(-a))``. The effective depth follows total material attenuation;
depth-independent pigment allocation and physiological yield cancel on conditioning.
This effective exponential closure is not an exact angular solution inside a leaf.
"""
struct ExponentialBirthProfile{FT<:AbstractFloat} <: AbstractLeafBirthProfile
    excitation_depth::FT
    forward_fraction::FT
    backward_fraction::FT
    function ExponentialBirthProfile{FT}(
        a::FT,
        forward::FT,
        backward::FT,
    ) where {FT<:AbstractFloat}
        isfinite(a) && a >= zero(FT) ||
            throw(ArgumentError("excitation depth must be finite and nonnegative"))
        isfinite(forward) && isfinite(backward) && min(forward, backward) >= zero(FT) ||
            throw(ArgumentError("birth amplitudes must be finite and nonnegative"))
        scale = max(forward, backward)
        scale > zero(FT) ||
            throw(ArgumentError("at least one birth amplitude must be positive"))
        f, b = forward / scale, backward / scale
        total = f + b
        new{FT}(a, f / total, b / total)
    end
end
ExponentialBirthProfile(a::FT, f::FT, b::FT) where {FT<:AbstractFloat} =
    ExponentialBirthProfile{FT}(a, f, b)

"""
    RegularPlateEscape(emission_depth, interface_transmission)

Symmetric regular plate with finite effective emission optical depth `b ≥ 0` and
the internal-interface-to-outside photon transmission probability `0 < t ≤ 1`.
The within-plate traversal survival is `exp(-b)`; it is distinct from `t`.
Both boundaries have reflection `1-t`. Emission starts equally toward either
face. Reflection paths end at the front face, back face or material reabsorption.
Perfectly trapping interfaces (`t=0`) are excluded from this three-fate contract.

Preparation/validation occurs on the host. The scalar application is allocation
free and supports Float32/Float64 CPU/CUDA; no optical database, intrinsic emission
spectrum, physiological yield, grid or device is owned by this object.
"""
struct RegularPlateEscape{FT<:AbstractFloat}
    emission_depth::FT
    interface_transmission::FT
    function RegularPlateEscape{FT}(b::FT, t::FT) where {FT<:AbstractFloat}
        isfinite(b) && b >= zero(FT) ||
            throw(ArgumentError("emission depth must be finite and nonnegative"))
        isfinite(t) && zero(FT) < t <= one(FT) ||
            throw(ArgumentError("interface transmission must be finite and in (0,1]"))
        new{FT}(b, t)
    end
end
RegularPlateEscape(b::FT, t::FT) where {FT<:AbstractFloat} = RegularPlateEscape{FT}(b, t)

# Integral of exp(-x*z) on the unit interval, including the removable singularity.
@inline emission_depth_integral(x) = iszero(x) ? one(x) : -expm1(-x) / x

# Normalized overlap integrals S(a,b)/L(a), D(a,b)/L(a). Scale sums before
# division, and combine exponentially small factors before taking exp. This
# also handles finite a,b whose sum overflows or whose separate L's underflow.
@inline function emission_overlaps(a::FT, b::FT) where {FT}
    scale = max(a, b)
    iszero(scale) && return (one(FT), one(FT))
    scaled_sum = a / scale + b / scale
    sum_loss = -expm1(-(scale * scaled_sum))
    if a > one(FT)
        a_loss = -expm1(-a)
        same = (a / scale) / scaled_sum * (sum_loss / a_loss)
        difference = abs(a - b)
        near = min(a, b)
        opposite = if difference > one(FT)
            exp(log(a / difference) - near) * (-expm1(-difference) / a_loss)
        else
            exp(log(a) - near) * (emission_depth_integral(difference) / a_loss)
        end
    else
        norm = emission_depth_integral(a)
        same =
            b <= one(FT) ? emission_depth_integral(a + b) / norm :
            (sum_loss / scaled_sum) / scale / norm
        opposite = exp(-min(a, b)) * emission_depth_integral(abs(a - b)) / norm
    end
    return same, opposite
end

# First four moments of the normalized forward excitation profile. Positive
# series avoids cancellation near zero; recurrence is well-conditioned at a>8.
@inline function emission_birth_moments(a::FT) where {FT}
    boundary = exp(-a) / emission_depth_integral(a)
    if a > FT(8)
        m1 = (one(FT) - boundary) / a
        m2 = (FT(2) * m1 - boundary) / a
        m3 = (FT(3) * m2 - boundary) / a
        m4 = (FT(4) * m3 - boundary) / a
        return m1, m2, m3, m4
    end
    term = (inv(FT(2)), inv(FT(3)), inv(FT(4)), inv(FT(5)))
    total = term
    for j = 1:64
        term = (
            term[1] * a / FT(j + 2),
            term[2] * a / FT(j + 3),
            term[3] * a / FT(j + 4),
            term[4] * a / FT(j + 5),
        )
        total = total .+ term
        maximum(term) <= eps(FT) * minimum(total) && break
    end
    return total .* boundary
end

# First-flight loss is symmetric under z -> 1-z, so it is independent of the
# forward/back mixture. Its weak-b series retains tiny true absorption instead
# of subtracting nearly unit escape probabilities. Constants are Taylor terms.
@inline function emission_first_loss(a::FT, b::FT, survival::FT) where {FT}
    if b < sqrt(sqrt(eps(FT)))
        m1, m2, m3, m4 = emission_birth_moments(a)
        s2 = one(FT) - FT(2) * m1 + FT(2) * m2
        s3 = one(FT) - FT(3) * m1 + FT(3) * m2
        s4 = one(FT) - FT(4) * m1 + FT(6) * m2 - FT(4) * m3 + FT(2) * m4
        return b * (inv(FT(2)) + b * (-s2 / FT(4) + b * (s3 / FT(12) - b * s4 / FT(48))))
    end
    return one(FT) - survival
end

"""
    leaf_emission_fates(escape, birth)

Return `(front, back, reabsorbed)` probabilities per already-produced internal
photon. Birth orientation uses fixed anatomical faces. Reabsorption is terminal;
no secondary emission, wavelength redistribution, absorptance or yield is applied.
Photon conservation holds to floating-point roundoff without clipping or rescaling.
Hemispheric face totals do not specify directional canopy radiance.

Built-in methods require the escape and exponential birth profile to share FT;
`UniformBirthProfile()` has no numeric parameters. Additional escape or birth
models can define methods without altering a LeafModel or host grid.
"""
function leaf_emission_fates end

@inline function regular_plate_fates(
    escape::RegularPlateEscape{FT},
    a::FT,
    f::FT,
    g::FT,
) where {FT}
    b, transmission = escape.emission_depth, escape.interface_transmission
    half = inv(FT(2))
    iszero(b) && return (front = half, back = half, reabsorbed = zero(FT))
    same, opposite = emission_overlaps(a, b)
    first_front = (f * same + g * opposite) * half
    first_back = (g * same + f * opposite) * half
    survival = (same + opposite) * half
    first_loss = emission_first_loss(a, b, survival)
    reflection = one(FT) - transmission
    loss = -expm1(-b)
    reflected_traversal = reflection * exp(-b)
    absorption_path = reflection * loss
    termination = transmission + absorption_path
    # Form the conditional branch before multiplying: t may be subnormal.
    escape_branch = transmission / termination
    front =
        (first_front + reflected_traversal * first_back) *
        (escape_branch / (one(FT) + reflected_traversal))
    back =
        (first_back + reflected_traversal * first_front) *
        (escape_branch / (one(FT) + reflected_traversal))
    reabsorbed = first_loss + survival * (absorption_path / termination)
    return (; front, back, reabsorbed)
end

@inline leaf_emission_fates(
    escape::RegularPlateEscape{FT},
    ::UniformBirthProfile,
) where {FT} = regular_plate_fates(escape, zero(FT), one(FT), zero(FT))
@inline leaf_emission_fates(
    escape::RegularPlateEscape{FT},
    birth::ExponentialBirthProfile{FT},
) where {FT} = regular_plate_fates(
    escape,
    birth.excitation_depth,
    birth.forward_fraction,
    birth.backward_fraction,
)
