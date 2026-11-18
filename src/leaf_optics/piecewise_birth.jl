@doc raw"""
    PiecewiseUniformBirthProfile(weights)

Birth density that is uniform within each of `N` equal-width layers of normalized
depth ``0 \le z \le 1`` (``z = x/d`` in a homogeneous plate of thickness ``d``,
``z = 0`` at the fixed anatomical front face, the coordinate of
[`ExponentialBirthProfile`](@ref)). Layer ``k`` spans
``[(k-1)/N,\, k/N]`` and receives the fraction ``w_k`` of the born photons:

```math
p(z) = N\,w_k \quad \text{for } (k-1)/N < z < k/N, \qquad \sum_{k=1}^{N} w_k = 1 .
```

`weights` is an `NTuple{N,FT}` or a host `AbstractVector{FT}` of finite nonnegative
relative weights with a positive sum; `N ≥ 1`. The host constructor stores them
normalized (scaled by their maximum, then divided by their compensated sum, so large
finite weights do not overflow and the stored fractions sum to one within about two
roundings for any `N`) in the field `layer_fractions::NTuple{N,FT}`; the object is
`isbits`.
Zero weights are allowed and describe dark layers. Invalid input throws
`ArgumentError`; mixed or non-floating-point element types are a `MethodError`.

This type carries a state-dependent birth profile prepared by the caller, typically
excitation fractions from [`birth_layer_fractions`](@ref) multiplied by the caller's
own layer allocation and yield (photons produced per layer). It owns no allocation,
yield, rate constant or excitation model. Because the plate is homogeneous, the
excitation and emission optical depths to depth ``z`` are ``a z`` and ``b z`` at
every wavelength; equal-width layers with uniform production inside each are a
discretization of the source, not tissue layers with their own optical properties.
The method for [`RegularPlateEscape`](@ref) is closed form; there is no method for
the two-zone [`DorsiventralPlateEscape`](@ref).
"""
struct PiecewiseUniformBirthProfile{N,FT<:AbstractFloat} <: AbstractLeafBirthProfile
    layer_fractions::NTuple{N,FT}
    function PiecewiseUniformBirthProfile{N,FT}(
        weights::NTuple{N,FT},
    ) where {N,FT<:AbstractFloat}
        N >= 1 || throw(ArgumentError("a piecewise birth profile needs at least one layer"))
        all(w -> isfinite(w) && w >= zero(FT), weights) ||
            throw(ArgumentError("layer weights must be finite and nonnegative"))
        scale = maximum(weights)
        scale > zero(FT) || throw(ArgumentError("at least one layer weight must be positive"))
        # Compensated sum: the stored fractions sum to one within about two
        # roundings for any N (a sequential sum would drift by about N roundings).
        total = compensated_sum(k -> weights[k] / scale, Val(N))
        new{N,FT}(ntuple(k -> (weights[k] / scale) / total, Val(N)))
    end
end
# A nonempty tuple binds FT; the empty tuple is rejected by its own method.
PiecewiseUniformBirthProfile(weights::Tuple{FT,Vararg{FT,M}}) where {FT<:AbstractFloat,M} =
    PiecewiseUniformBirthProfile{M + 1,FT}(weights)
PiecewiseUniformBirthProfile(::Tuple{}) =
    throw(ArgumentError("a piecewise birth profile needs at least one layer"))
function PiecewiseUniformBirthProfile(weights::AbstractVector{FT}) where {FT<:AbstractFloat}
    isempty(weights) &&
        throw(ArgumentError("a piecewise birth profile needs at least one layer"))
    return PiecewiseUniformBirthProfile(Tuple(weights))
end

# Neumaier's compensated summation step: the rounding error of every addition is
# recovered exactly and carried in `compensation`. No branch on data beyond a
# select, so it is GPU friendly; the error bound of the final sum + compensation is
# about two roundings independently of the number of terms.
@inline function compensated_add(total::FT, compensation::FT, term::FT) where {FT}
    sum = total + term
    compensation += abs(total) >= abs(term) ? (total - sum) + term : (term - sum) + total
    return sum, compensation
end

@inline function compensated_sum(term, ::Val{N}) where {N}
    total, compensation = compensated_add(zero(term(1)), zero(term(1)), term(1))
    for k = 2:N
        total, compensation = compensated_add(total, compensation, term(k))
    end
    return total + compensation
end

# 1 - L(x) = (x + expm1(-x))/x without cancellation. Below one, the alternating
# series sum_{n>=1} (-1)^(n+1) x^n/(n+1)! (terms falling by x/(n+2), at most 19
# terms in Float64). From one on, x + expm1(-x) >= exp(-1), so the sum loses at
# most a factor x/(x + expm1(-x)) <= e to cancellation.
@inline function emission_depth_deficit(x::FT) where {FT}
    if x < one(FT)
        term = x / FT(2)
        total = term
        for n = 2:40
            term *= -x / FT(n + 1)
            total += term
            abs(term) <= eps(FT) / 2 * total && break
        end
        return total
    end
    return (x + expm1(-x)) / x
end

@doc raw"""
    leaf_emission_fates(escape::RegularPlateEscape{FT}, birth::PiecewiseUniformBirthProfile{N,FT})

Closed-form fates of a regular plate for a layer-wise uniform birth profile. With
``x = b/N``, ``L(x) = -\operatorname{expm1}(-x)/x`` and layer weights ``w_k``, the
first-flight probabilities of reaching the front and back interface are

```math
F_1 = \tfrac12 L(x)\sum_{k=1}^{N} w_k\, e^{-b (k-1)/N}, \qquad
B_1 = \tfrac12 L(x)\sum_{k=1}^{N} w_k\, e^{-b (N-k)/N},
```

and the first-flight reabsorption is formed from nonnegative terms only,

```math
A_1 = \tfrac12 \sum_{k=1}^{N} w_k \Bigl[m_{k-1} + m_{N-k}
      + \bigl(e^{-b(k-1)/N} + e^{-b(N-k)/N}\bigr)\bigl(1 - L(x)\bigr)\Bigr],
\qquad m_j = -\operatorname{expm1}(-b j/N),
```

with ``1 - L(x)`` from a series below ``x = 1``, and ``A_1 = b/2`` for
``b \le \varepsilon_{\mathrm{FT}}`` (exact to relative ``O(b)`` for any profile). The three layer sums use compensated
(Neumaier) summation, so their accuracy does not degrade with `N`. The
multiple-reflection closure is the one of the exponential profile. Cost is one
exponential per layer.
"""
@inline function leaf_emission_fates(
    escape::RegularPlateEscape{FT},
    birth::PiecewiseUniformBirthProfile{N,FT},
) where {N,FT}
    b = escape.emission_depth
    half = inv(FT(2))
    iszero(b) && return (front = half, back = half, reabsorbed = zero(FT))
    w = birth.layer_fractions
    x = b / FT(N)
    mean_survival = emission_depth_integral(x)
    deficit = emission_depth_deficit(x)
    # Compensated accumulators (sum, compensation): with nonnegative terms the
    # first-flight sums are accurate to a few roundings for any N.
    toward_front = toward_front_c = zero(FT)
    toward_back = toward_back_c = zero(FT)
    boundary_loss = boundary_loss_c = zero(FT)
    # One transcendental per boundary: the complement of the smaller of exp(-c)
    # and -expm1(-c) is formed by subtraction without cancellation.
    half_survival_depth = log(FT(2))
    # Index j is a layer boundary at depth j/N from one face: the near edge of layer
    # j+1 seen from the front and of layer N-j seen from the back.
    for j = 0:N-1
        c = x * FT(j)
        if c <= half_survival_depth
            lost = -expm1(-c)
            transmitted = one(FT) - lost
        else
            transmitted = exp(-c)
            lost = one(FT) - transmitted
        end
        front_weight, back_weight = w[j+1], w[N-j]
        toward_front, toward_front_c =
            compensated_add(toward_front, toward_front_c, front_weight * transmitted)
        toward_back, toward_back_c =
            compensated_add(toward_back, toward_back_c, back_weight * transmitted)
        # (front + back weight) keeps the loss symmetric under reversed weights.
        boundary_loss, boundary_loss_c =
            compensated_add(boundary_loss, boundary_loss_c, (front_weight + back_weight) * lost)
    end
    toward_front += toward_front_c
    toward_back += toward_back_c
    boundary_loss += boundary_loss_c
    first_front = half * mean_survival * toward_front
    first_back = half * mean_survival * toward_back
    survival = first_front + first_back
    # For b <= eps the first-flight loss is b/2 to relative O(b) for any birth
    # profile (half of b z + b (1 - z)). This also avoids x = b/N rounding on the
    # subnormal grid for b < N floatmin, where the layer form loses about N/4
    # subnormal units.
    first_loss = b <= eps(FT) ? half * b :
                 half * (boundary_loss + deficit * (toward_front + toward_back))
    return regular_plate_reflections(escape, first_front, first_back, survival, first_loss)
end

@doc raw"""
    birth_layer_fractions(birth::ExponentialBirthProfile{FT}, Val(N)) -> NTuple{N,FT}
    birth_layer_fractions(birth, Val(N), FT) -> NTuple{N,FT}

Fractions of the normalized birth (excitation) density in `N` equal-width layers of
dimensionless depth, the layer coordinate of [`PiecewiseUniformBirthProfile`](@ref).
For [`ExponentialBirthProfile`](@ref) with depth ``a`` and normalized amplitudes
``(f, g)`` the exact layer integrals are

```math
P_k = \frac{\operatorname{expm1}(-a/N)}{\operatorname{expm1}(-a)}
      \Bigl[f\,e^{-a(k-1)/N} + g\,e^{-a(N-k)/N}\Bigr],
```

evaluated as ``L(a/N)/(N L(a))`` times the bracket for ``a < 1`` so that ``a \to 0``
gives ``1/N``, and as the ratio of `expm1` otherwise so that large finite ``a`` neither
overflows nor underflows the prefactor. They sum to one to roundoff.

The three-argument form takes the floating-point type explicitly. It is required for
[`UniformBirthProfile`](@ref), which carries no number and gives `N` copies of
`1/N`, and is accepted for an exponential profile of the same `FT` (another `FT` is a
`MethodError`), so generic code can always pass it. `N ≥ 1` is a compile-time
constant; the result is `isbits` and the evaluation allocation free.

This is the physics-free half of a depth-dependent source: the caller multiplies
the fractions by its own antenna excitation, layer allocation to the source and
layer yield, sums the products for the source's production and passes the products
to [`PiecewiseUniformBirthProfile`](@ref). No allocation or yield is supplied here.
"""
@inline function birth_layer_fractions(
    birth::ExponentialBirthProfile{FT},
    ::Val{N},
) where {FT,N}
    N isa Int && N >= 1 ||
        throw(ArgumentError("number of layers must be a positive integer"))
    a, f, g = birth.excitation_depth, birth.forward_fraction, birth.backward_fraction
    x = a / FT(N)
    share = if a < one(FT)
        emission_depth_integral(x) / (FT(N) * emission_depth_integral(a))
    else
        expm1(-x) / expm1(-a)
    end
    return ntuple(k -> share * (f * exp(-x * FT(k - 1)) + g * exp(-x * FT(N - k))), Val(N))
end
@inline birth_layer_fractions(
    birth::ExponentialBirthProfile{FT},
    layers::Val,
    ::Type{FT},
) where {FT<:AbstractFloat} = birth_layer_fractions(birth, layers)
@inline function birth_layer_fractions(
    ::UniformBirthProfile,
    ::Val{N},
    ::Type{FT},
) where {N,FT<:AbstractFloat}
    N isa Int && N >= 1 ||
        throw(ArgumentError("number of layers must be a positive integer"))
    return ntuple(_ -> inv(FT(N)), Val(N))
end
