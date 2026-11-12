# Angular normalization for scattering and emission

CanopyOptics evaluates a continuous leaf-angle distribution with quadrature.
That is distinct from the RT solver's quadrature over outgoing light directions.
A continuously normalized Lambertian leaf law need not integrate to exactly one
on a finite outgoing-angle rule.

The default remains [`ContinuousAngularNormalization`](@ref). The opt-in
[`DiscreteAngularNormalization`](@ref) makes each leaf face's redistribution
conservative on the supplied angular measure before integrating over leaf
inclination. It leaves leaf reflectance/transmittance, G and leaf-area weights
unchanged. It does not compute fluorescence yield, change APAR after a solve,
or supply a new transport solver.

## Use the same measure for scattering and a leaf source

The RT backend supplies positive-hemisphere direction cosines `μ` and weights
`weights` integrating `dμ` over `[0,1]`. These are **not** the flux weights `μ*w`
and must not be repeated across Stokes components. Observation/solar directions
may have zero weight. Numerical tolerances and physical parameters are explicit
caller inputs:

```julia
normalization = DiscreteAngularNormalization(μ, weights; atol=angular_tolerance)
leaf = BiLambertianCanopyScattering(R, T)
quadrature = CanopyQuadrature(n_leaf=leaf_nodes)

Zpp, Zmp = compute_Z_matrices(leaf, μ, leaf_distribution, 0:m_max;
    quadrature, angular_normalization=normalization)

# Reuse this same numerical measure when constructing escaped leaf emission.
P, N = leaf_projection_moments(ray_cosine, leaf_normal_cosine, m_max;
    normalization)
source_moments = 2 .* (P .* front_escape .+ N .* back_escape)
```

The last expression precedes the leaf-normal azimuth factor and any adapter
conversion to the solver's azimuth convention. The solver's source path factor
`ΔLAI/abs(μ)` is also separate. `front_escape` and `back_escape` already represent
photons that escaped the leaf; do not apply leaf absorption or yield again.

[`leaf_angular_scale`](@ref) exposes the host geometry factor for adapters that
already cache their projection matrices. Apply it once to each leaf orientation,
for both faces and every Fourier mode. Build geometry on the host, then copy
optical/source arrays to the desired CPU/GPU backend. This API does not add a GPU
RT orchestration layer.

The constructor owns copies of the angular arrays. Treat them as read-only.
Scattering assembly rejects mismatched angle values, ordering and precision.
Rebuild the normalization and downstream optical/source caches whenever angles
or weights change. A cache key must include the normalization choice and the
values of both angular arrays; matching dimensions alone is insufficient.

## Conservation and reciprocity

Let `P₀` and `N₀` be the positive/negative zeroth projection moments of a leaf
whose inclination cosine is `μ_leaf`. The discrete integral of a leaf face is

```math
H_l = 2\sum_i w_i(P_{0,il}+N_{0,il}).
```

The continuous value is one. The discrete convention multiplies that leaf's
scattering contribution by `1/H_l`. Because one scalar weights the entire
incoming/outgoing projection product for that leaf, the scalar area kernel
`Z * Diagonal(G)` remains reciprocal. The zeroth-mode outgoing integral becomes

```math
\frac12\sum_iw_i(Z^{++}_{ij}+Z^{-+}_{ij})=1
```

up to floating-point errors. This also normalizes each front/back emitted source
on the same discrete measure. It is different from dividing assembled columns
by their individual errors, which need not preserve reciprocity, or rescaling
a completed flux/APAR calculation to hide a budget residual.

Finite angular approximation still matters. On a symmetric Gaussian grid the
front-face upward fraction changes from the continuous `(1+μ_leaf)/2` to
`(1+μ_leaf/H_l)/2`; this difference converges away with angular refinement.
Check viewing radiances and hemispheric partitioning, not only conservation.
Finite elemental-depth and roundoff errors also remain in a downstream solver.

## Scope and validation

The option is implemented for the bi-Lambertian analytic, single-moment and
range APIs, plus their compatibility wrappers. Its `npol` expansion remains
scalar diffuse I→I; this is not a polarized/specular conservation treatment.
Other canopy-component APIs retain their existing behavior and do not silently
apply this policy. Both Float32 and Float64 are tested. Existing differentiation
with respect to leaf optical parameters is retained; differentiation of the
discrete angular-grid construction is not promised.

Tests cover column conservation, reciprocity across Fourier modes, exact
horizontal/vertical leaf limits, source photon budgets, zero-weight viewing
nodes, grid mismatch rejection, owned input arrays, black/reflecting/transmitting
limits, four continuous beta leaf-angle PDFs, Stokes expansion and optical
parameter differentiation. Default calls and explicit continuous normalization
are checked for equality. A separate CUDA acceptance script checks host-built
source geometry transferred to device arrays.
