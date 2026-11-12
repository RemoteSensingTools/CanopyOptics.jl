# Scalar plate algebra. No clipping or epsilon absorption is introduced.
# `prospect_primal` is defined in `src/types/canopy_types.jl`.

"""Internal transmission and absorptance for diffuse illumination of one slab."""
@inline function prospect_internal_attenuation(k)
    value = prospect_primal(k)
    isfinite(value) && value >= zero(value) ||
        throw(ArgumentError("Material optical depth must be finite and nonnegative"))
    if value <= one(value)
        # 1-2E3(k), expressed without subtracting almost equal numbers.
        # At zero this retains a'(0)=2 and tau'(0)=-2, including Dual inputs.
        # The physical second derivative is singular at k=0; no Hessian is claimed.
        a = iszero(value) ? k+k : -expm1(-k) + k*exp(-k) - k*k*expint(k)
        return (; transmission = one(k)-a, absorptance = a)
    end
    # Independent transmission avoids losing a small positive tau in 1-a.
    # Keep the generalized exponential-integral order an integer for AD.
    # SpecialFunctions' general-order implementation has mixed FT return paths.
    # Convert its result (including Dual tangents) back to the caller's type.
    tau = (one(k)+one(k))*oftype(k, expint(3, k))
    return (; transmission = tau, absorptance = one(k)-tau)
end

"""
Fractional stack of m>=0 identical diffuse plates. `alpha` is computed local
absorption, not the residual of rounded reflection/transmission. The analytic
near-conservative series preserves first derivatives at alpha=0 and m=0.
"""
@inline function prospect_plate_stack(r, t, alpha, m)
    if iszero(prospect_primal(t))
        # Opaque-layer and zero-thickness limits do not commute. At their joint
        # endpoint there is no finite general derivative with respect to m.
        if iszero(prospect_primal(m))
            return (; transmission = one(t), reflectance = zero(t), absorptance = zero(t))
        end
        return (; transmission = zero(t), reflectance = r, absorptance = alpha)
    end
    two=one(t)+one(t)
    h=alpha*(one(t)-t+r)
    # Decide the series domain without first forming h/t, which can overflow
    # for a positive subnormal transmission. Numerical threshold, not a loss floor.
    scale=max(one(prospect_primal(m)), abs(prospect_primal(m)))
    threshold=cbrt(eps(typeof(prospect_primal(t))))
    small=iszero(prospect_primal(h)) ||
          prospect_primal(h)*scale <= two*prospect_primal(t)*(threshold/scale)
    if small
        z=h/(two*t)
        w=(m*(h/t))*(m/two)
        # Chebyshev continuation C=cosh(m acosh(1+z)), S=dC/dz / m.
        # Through z^3 in C-1 and z^3 in S; omitted terms are O((m^2 z)^4).
        cminus=w*(one(z)+(w-z)/6+(w-z)*(w-4z)/90)
        s=m*(one(z) + (w-z)/3 + (w-z)*(w-4z)/30 + (w-z)*(w-4z)*(w-9z)/630)
        denominator=t*(one(z)+cminus)+(one(t)-t-h/two)*s
        return (;
            transmission = t/denominator,
            reflectance = r*s/denominator,
            absorptance = (t*cminus+alpha*(one(t)+t-r)*s/two)/denominator,
        )
    end
    d=sqrt(h)*sqrt(h+4t)
    mu=prospect_primal(h)<=prospect_primal(t) ? log1p((h/two+d/two)/t) :
       log(t+h/two+d/two)-log(t)
    argument=m*mu
    q=exp(-argument)
    u=-expm1(-two*argument)
    c=(one(q)+q*q)/two
    cminus=expm1(-argument)^2/two
    ratio=u/d
    boundary=one(t)-t-h/two
    denominator=c+boundary*ratio
    # For thick/weakly absorbing stacks, the quotient r*ratio/denominator
    # cancels large terms in its derivative. Divide numerator and denominator
    # by ratio first; retain the original form near m=0 to preserve that tangent.
    reflectance=prospect_primal(ratio)>one(prospect_primal(ratio)) ? r/(c/ratio+boundary) :
                r*ratio/denominator
    return (;
        transmission = q/denominator,
        reflectance,
        absorptance = (cminus+alpha*(one(t)+t-r)*ratio/two)/denominator,
    )
end

"""One spectral PROSPECT plate system with independently accumulated absorption."""
@inline function prospect_spectral_budget(k, N, talf, t12, nr)
    internal=prospect_internal_attenuation(k)
    tau, a=internal.transmission, internal.absorptance
    t21=t12/(nr*nr)
    r21=one(t21)-t21
    boundary=t21+r21*a
    denominator=boundary*(one(tau)+r21*tau)
    t=t12*tau*t21/denominator
    alpha=t12*a/boundary
    # Reflection follows the same independent local absorption/transmission
    # solution; no absorbing result is rescaled to force a flux budget.
    r=one(t)-t-alpha
    ta=talf*tau*t21/denominator
    aa=talf*a/boundary
    ra=one(ta)-ta-aa
    remainder=prospect_plate_stack(r, t, alpha, N-one(N))
    coupling=one(r)-remainder.reflectance*r
    transmission=ta*remainder.transmission/coupling
    reflectance=ra+ta*remainder.reflectance*t/coupling
    absorptance=aa+ta*(remainder.absorptance+remainder.reflectance*alpha)/coupling
    return (; transmission, reflectance, absorptance)
end
