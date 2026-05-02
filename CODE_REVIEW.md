# CanopyOptics.jl — Code Review

Review of uncommitted changes on `main` as of 2026-05-01.
Reviewed files: `Project.toml`, `src/CanopyOptics.jl`,
`src/core/leafAngleRoutines.jl`, `test/test_CanopyOptics.jl`.

---

## Summary of changes

The diff has three logically distinct pieces:

1. **Bug fix.** Both overloads of `compute_Z_matrices_aniso` previously wrote
   `Z[i_in, :]` (rows = incoming) and omitted the `1/G(μ)` normalisation.
   Both now write `Z[:, i_in]` and divide by `G[i] · ϖ`, matching the
   canonical canopy convention (`Z[i_out, j_in]`, with rows = outgoing μ).
2. **New analytic kernel.** `compute_Z_matrices_aniso_analytic` replaces the
   brute-force outgoing-azimuth quadrature with a closed-form Fourier
   expansion of the clipped leaf-projection function
   (`_clipped_projection_moments`). Faster, no azimuth quadrature noise,
   produces all `m = 0:m_max` moments in one call.
3. **Type-stability cleanups in `G`** plus removal of stray `@show` debug
   statements.

The math in (2) is sound and the new test set in `test_CanopyOptics.jl`
covers it well: P₀ ≡ Shultis–Myneni `H`, exact `P_m − N_m` identity,
agreement with the old brute-force code, flux-conservation column sum
≈ 2, reciprocity `Z·G = (Z·G)ᵀ`, smoke test up to m=64.

---

## Round 1 — Correctness / blocking issues

### 1. `Project.toml` compat — RESOLVED

Earlier draft of the diff had `LazyArtifacts = "1.11.0"` (which would have
broken Julia ≤ 1.10). Current diff uses `LazyArtifacts = "1"`, with the
package added to both `[deps]` and `[compat]`, and `julia` compat updated
to include `1.12`. ✅

### 2. `compute_Z_matrices_aniso(...,Zup,Zdown,m)` is a silent behaviour change

The precomputed-Z overload at `src/core/leafAngleRoutines.jl:562` used to
write `𝐙[i,:]`; now it writes `𝐙[:,i]` **and** divides by `G[i]`.
Anything outside this repo using the old behaviour will now get a
transposed, differently-scaled matrix without warning. The local consumer
is `vSmartMOM.jl/src/CoreRT/Surfaces/canopy_surface.jl` (also uncommitted),
so it must be updated in lockstep — and a one-line release note ("output
convention now `Z[i_out, j_in]` and divided by `G·ϖ`") would save a
future user a debugging session.

### 3. `G` test tolerance loosened from `0.001` → `0.002`

```julia
- @test all(abs.(G .- 0.5) .< 0.001)
+ @test all(abs.(G .- 0.5) .< 0.002)
```

Two possibilities, only one of which is acceptable:

- The old test was passing because some quantities were silently promoted
  to `Float64` via `2θₗ/π` (Float64 `π`), and the new strictly-FT-typed
  version loses ~1 ULP. Fine for `Float32`; for `Float64` you should
  still be at `< 0.001`.
- Or the change accidentally hurt accuracy.

Prefer:

```julia
@test all(abs.(G .- 0.5) .< (FT == Float32 ? 2e-3 : 1e-3))
```

…rather than blanket-loosening for both float types.

### 4. Repeated allocation in the analytic hot loop

```julia
for l in eachindex(θₗ)
    for i in eachindex(μ)
        P, N   = _clipped_projection_moments( μ[i], μ_L, m_max)
        Pm, Nm = _clipped_projection_moments(-μ[i], μ_L, m_max)
        ...
```

Each call allocates two `Vector{FT}` of length `m_max + 1`. Total
allocations ≈ `4 · nQuad · nμ` vectors per assembly. For typical sizes
(nQuad=64, nμ=8, m_max=33) that's ~2000 small allocations. Removable by
passing pre-allocated buffers — see Round 2 §5 (mutating-variant
pattern).

### 5. Loop traversal order

```julia
for k in 1:nm, j in 1:nμ, i in 1:nμ
    @inbounds begin
        Ψpp_same = _projection_same_vsmartmom(Pꜜ[j, k], Nꜜ[j, k], Pꜜ[i, k], Nꜜ[i, k])
        ...
        Z⁺⁺[i, j, k] += leaf_weight * (...)
```

Innermost is `i`, indexing `Z⁺⁺[i, j, k]` — column-major friendly, good.
But `Pꜜ[j, k]` is loop-invariant inside the `i` loop; the compiler
probably hoists it, but explicit hoisting (`Pdj = Pꜜ[j, k]` before the
`i` loop) is one less hope.

### 6. Docstring nitpicks

- The `_clipped_projection_moments` docstring has `function `H``
  (single-quote-then-backtick-`H`) which renders awkwardly and pollutes
  `grep -n "^function"` output. Cosmetic.
- Module-level convention block in `CanopyOptics.jl` is excellent — keep.
- The `compute_Z_matrices_aniso_analytic` docstring's WHY note about the
  leading `2` ("converts the `1/(2π)` Fourier coefficients to the
  half-range cosine moment used by the historical implementation") is
  exactly the right kind of comment. Keep.

### 7. Minor

- `_canopy_fourier_factor` is defined but only used by
  `_normalise_vsmartmom_Z!`. See Round 2 §6 — inline it.
- `compute_Z_matrices_aniso_analytic` requires `ϖ > 0`. Throwing on the
  fully-absorbing-leaf limit (`R = T = 0`) is unfriendly; return zero
  matrices instead. See Round 2 §10.
- `_leaf_inclination_quadrature` returns `(θₗ, w .* F)` already
  pre-multiplied. Future contributor will assume `wₗ` is just Gauss
  weights and apply F twice. See Round 2 §9.

---

## Round 1 — What looks right

- `_clipped_projection_moments` correctly handles all three branches
  (`a ≥ b`, `a ≤ −b`, `|a| < b`) and the degenerate `b == 0` case, and
  the `P − N = a δ_{m0} + (b/2) δ_{m1}` identity is exact (test
  verifies this to `1e-14`).
- `m=0` matches the existing `compute_Z_matrices` (Shultis–Myneni Eq. 45
  assembly via `compute_Ψ`) to `1e-12` — strong cross-check that the
  closed-form Fourier moments reduce correctly.
- Removing the `@show` calls in `bfG`, `compute_Γ`, `precompute_Zazi_`
  is a clean win.
- The reciprocity test (`Z·G == (Z·G)ᵀ`) is the right invariant to
  assert — structural property of the analytic formula, doesn't ride on
  quadrature accuracy.

### Suggested commit shape

Split into three commits for easier review/bisect:

1. `fix(Project.toml): correct LazyArtifacts compat and re-sort`
2. `fix(Z_aniso): use canonical convention (Z[:, i_in] / G·ϖ) in both overloads`
3. `feat(canopy): closed-form Fourier moments for BiLambertian Z (compute_Z_matrices_aniso_analytic)`

Tests land alongside (3).

---

## Round 2 — Julian style, dispatch, long-term maintainability

None of the items below are blocking; they're about whether this code
will age well.

### 1. Don't bake a downstream consumer's name into upstream library helpers

`_projection_same_vsmartmom`, `_projection_opposite_vsmartmom`,
`_normalise_vsmartmom_Z!` — three private helpers in `CanopyOptics.jl`
named after `vSmartMOM`. The convention they implement (`Z[i_out, j_in]`,
divide by `G·ϖ`, half-range cosine moments) is a *canopy* convention;
vSmartMOM happens to share it because that's what the canonical
Shultis–Myneni / Knyazikhin–Marshak algebra produces. If a year from now
someone uses `CanopyOptics` from a non-vSmartMOM RT solver, these names
will read as historical accidents.

**Suggested rename**:

| Current                          | Suggested        |
| -------------------------------- | ---------------- |
| `_projection_same_vsmartmom`     | `_psi_same`      |
| `_projection_opposite_vsmartmom` | `_psi_opposite`  |
| `_normalise_vsmartmom_Z!`        | `_normalise_Z!`  |

The big docstring on `compute_Z_matrices_aniso_analytic` already says
"vSmartMOM convention" — that's the right place to mention the consumer.

### 2. Two implementations of the same physics will drift

After this diff, the package has three routines computing canopy Z:

| Routine                                                     | What it does                                          |
| ----------------------------------------------------------- | ----------------------------------------------------- |
| `compute_Z_matrices(mod, μ, LD, m)`                         | m=0 only, via `compute_Ψ` (Shultis–Myneni Eq. 45)     |
| `compute_Z_matrices_aniso(mod, μ, LD, m)` (and `Zup,Zdown`) | brute-force outgoing-azimuth quadrature               |
| `compute_Z_matrices_aniso_analytic(mod, μ, LD, m_max)`      | closed-form Fourier moments                           |

Long term these will drift unless they share a backbone. Two refactors
worth the effort:

**(a) Make the old `compute_Z_matrices` (m=0) call the analytic path.**
The m=0 result of `compute_Z_matrices_aniso_analytic` *is* what
`compute_Z_matrices` computes — they should not be parallel
implementations of the same equation. Once verified equivalent (the
test already does this to 1e-12), redefine `compute_Z_matrices` as a
thin wrapper that slices the analytic output. Delete `compute_Ψ`'s
production usage.

**(b) Demote `compute_Z_matrices_aniso(mod, μ, LD, m)` (brute-force) to
a test reference.** It's the only place you exercise `compute_Γ`, the
per-direction integrand. After this diff, it has no production caller —
the analytic path supersedes it. Move it to `test/` (or keep it under a
`# Test reference` section header) so future readers don't mistake it
for the production API.

If both moves are done, the production surface collapses to one
function: `compute_Z_matrices_aniso_analytic`. That's the right Julian
outcome — one canonical method per physics question.

### 3. Multiple dispatch — the API name is doing the dispatching

`compute_Z_matrices` vs `compute_Z_matrices_aniso` vs
`compute_Z_matrices_aniso_analytic` — three names because the *return
shape* (2D vs 3D), the *backend* (brute-force vs analytic), and the
*isotropy assumption* (`compute_Z_matrices` ignores m≥1) all change with
the suffix. That's name-based dispatch, not type-based dispatch.

A more Julian factoring:

```julia
abstract type CanopyZBackend end
struct AnalyticFourier <: CanopyZBackend end
struct BruteForceAzimuth <: CanopyZBackend
    nQuad::Int = 64
end

# Single-m: returns (Z⁺⁺, Z⁻⁺) :: (Matrix, Matrix)
compute_Z(mod, μ, LD, m::Int; backend = AnalyticFourier()) = ...

# All-m up to m_max: returns (Z⁺⁺, Z⁻⁺) :: (3D, 3D)
compute_Z(mod, μ, LD, m_range::AbstractUnitRange; backend = AnalyticFourier()) = ...
```

This makes the backend choice explicit and the rank-of-output decision
driven by the `m` argument's type. You also get a natural place to add
a `GPUFourier` backend later without inventing a fourth function name.

If that's too much surgery, at minimum: the new function shouldn't be
named `compute_Z_matrices_aniso_analytic` (suffix soup). Either
`compute_Z_matrices_fourier` (signals what's different — closed-form
Fourier moments) or fold into `compute_Z_matrices` with a 4th positional
arg `m_max::Int`.

### 4. Type signature is more restrictive than it needs to be

```julia
function compute_Z_matrices_aniso_analytic(mod::BiLambertianCanopyScattering,
                                           μ::AbstractVector{FT},
                                           LD::AbstractLeafDistribution,
                                           m_max::Int) where {FT<:AbstractFloat}
```

Two concerns:

- **`AbstractFloat` blocks `ForwardDiff.Dual`.** The package adds
  `ForwardDiff = "0.10"` to `[compat]` and `prospect` already supports
  ForwardDiff. If someone wants `∂Z/∂R` or `∂Z/∂T`, they can't, because
  `Dual <: Real` but not `<: AbstractFloat`. Loosen to `<:Real` — the
  body uses `acos`, `sin`, `sqrt`, `clamp`, all of which dispatch
  through ForwardDiff fine.
- **`mod` and `μ` independently parametric.** Right now
  `BiLambertianCanopyScattering{FT}` forces `mod.R::FT, mod.T::FT`, but
  `μ` is also `FT` — they must match. If a user passes
  `BiLambertianCanopyScattering{Float32}` with `μ::Vector{Float64}`,
  this won't dispatch (silent `MethodError` rather than a useful
  message). Either:

  ```julia
  function compute_Z_matrices_aniso_analytic(mod::BiLambertianCanopyScattering,
                                             μ::AbstractVector{<:Real},
                                             LD::AbstractLeafDistribution,
                                             m_max::Int)
      FT = promote_type(eltype(μ), typeof(mod.R))
      ...
  ```

  …or document that they must match.

The same loosening applies to `_clipped_projection_moments` — which
currently has `where {FT<:Real}` (good!) but the
`compute_Z_matrices_aniso_analytic` caller is stricter than its callee.

### 5. Mutating-variant pattern for hot-path allocations

`_clipped_projection_moments` is the hottest function in the analytic
path; called `nQuad · nμ · 2` times per assembly. Julian convention: the
in-place version takes a `!` suffix and accepts buffers; the allocating
version is a thin wrapper.

```julia
function _clipped_projection_moments!(P::AbstractVector{FT}, N::AbstractVector{FT},
                                      μ::FT, μ_L::FT) where {FT<:Real}
    @assert length(P) == length(N)
    m_max = length(P) - 1
    fill!(P, zero(FT))
    fill!(N, zero(FT))
    # ... existing body, writing into P and N in place ...
    return P, N
end

function _clipped_projection_moments(μ::FT, μ_L::FT, m_max::Int) where {FT<:Real}
    P = Vector{FT}(undef, m_max + 1)
    N = Vector{FT}(undef, m_max + 1)
    _clipped_projection_moments!(P, N, μ, μ_L)
end
```

Then `compute_Z_matrices_aniso_analytic` allocates four scratch buffers
once outside the leaf-quadrature loop and reuses them. Standard Julia
stdlib pattern (`mul!`, `lmul!`, `LinearAlgebra` everywhere).

### 6. `_canopy_fourier_factor` is over-abstracted

```julia
@inline function _canopy_fourier_factor(::Type{FT}, m::Int) where {FT}
    return m == 0 ? FT(2) : FT(4)
end
```

Single-use, three-line function for `m == 0 ? 2 : 4`. Inline it into
`_normalise_vsmartmom_Z!`:

```julia
ff = (k == 1) ? FT(2) : FT(4)   # m = k-1
```

Premature abstraction; reserve named helpers for things that have a
*physics* name (`_psi_same`, `_clipped_projection_moments`).

### 7. Hot-loop kernel deserves its own function for testability and future GPU

The triple loop in `compute_Z_matrices_aniso_analytic`:

```julia
for k in 1:nm, j in 1:nμ, i in 1:nμ
    @inbounds begin
        Ψpp_same = _projection_same_vsmartmom(Pꜜ[j, k], Nꜜ[j, k], Pꜜ[i, k], Nꜜ[i, k])
        ...
        Z⁺⁺[i, j, k] += leaf_weight * (T_leaf * Ψpp_same + R_leaf * Ψpp_opp)
        Z⁻⁺[i, j, k] += leaf_weight * (T_leaf * Ψmp_same + R_leaf * Ψmp_opp)
    end
end
```

…is a textbook KernelAbstractions kernel: per-element work, no inner
reduction, all reads from arrays sliced by outer index `k`. Pulling it
out as a function:

```julia
function _accumulate_Z!(Z⁺⁺, Z⁻⁺, Pꜜ, Nꜜ, Pꜛ, Nꜛ, R, T, weight)
    nμ, nm = size(Pꜜ)
    for k in 1:nm, j in 1:nμ
        Pdj, Ndj = Pꜜ[j, k], Nꜜ[j, k]
        @inbounds for i in 1:nμ
            Pdi, Ndi = Pꜜ[i, k], Nꜜ[i, k]
            Pui, Nui = Pꜛ[i, k], Nꜛ[i, k]
            Z⁺⁺[i, j, k] += weight * (T * 2*(Pdj*Pdi + Ndj*Ndi) + R * 2*(Pdj*Ndi + Ndj*Pdi))
            Z⁻⁺[i, j, k] += weight * (T * 2*(Pdj*Pui + Ndj*Nui) + R * 2*(Pdj*Nui + Ndj*Pui))
        end
    end
end
```

…buys you (a) hoisted `Pdj, Ndj` loads (compiler may already do this,
but explicit is better), (b) testability on a single (j, μ_L), (c) a
future `@kernel` rewrite that's a straight transcription. The current
body interleaves the per-leaf-inclination quadrature with the per-(i,j)
accumulation, which makes the kernel boundary fuzzy.

### 8. Comment hygiene — historical bug commentary will rot

```julia
# NB: previous version of this function wrote rows-as-incoming
# (𝐙[i,:]) which produced a transposed result vs the single-shot
# path; now both write columns-as-incoming (𝐙[:,i]).
```

Fine commit-message line. As an in-source comment it's a pure historical
note that becomes confusing once "previous" is several years back. After
commit, leave only the *positive* convention statement:

```julia
# Convention: Z[i_out, j_in], rows = outgoing μ, columns = incoming μ.
```

### 9. `_leaf_inclination_quadrature` returns a measure-weighted weight

```julia
return θₗ, wθ .* Fₗ
```

The returned `wₗ` is `Gauss_weight × leaf_pdf × scaling`. A future
contributor will see `wₗ` at the call site and assume Gauss weights,
then multiply by `pdf(LD, ...)` again. Either rename:

```julia
return θₗ, wθ .* Fₗ   # second return = "leaf-measure quadrature weight" (incl. pdf)
```

…with that comment on the return line, or change the function's return
to a `NamedTuple`:

```julia
return (θ = θₗ, w_measure = wθ .* Fₗ)
```

…so the call site reads
`(θ, w_measure) = _leaf_inclination_quadrature(...)` and the surprise is
gone.

### 10. `compute_Z_matrices_aniso_analytic` argument check throws on absorbing leaf

```julia
ϖ <= zero(FT) && throw(ArgumentError("R + T must be positive ..."))
```

`R = T = 0` is a legitimate physical limit (perfectly absorbing leaves).
Throwing is unfriendly. Better: short-circuit return zero matrices:

```julia
ϖ = R_leaf + T_leaf
if ϖ ≤ zero(FT)
    return zeros(FT, nμ, nμ, nm), zeros(FT, nμ, nμ, nm)
end
```

The downstream RT solver will multiply by ϖ anyway, so zeros are
correct, not an error.

### 11. Symmetry: `_clipped_projection_moments(-μ, μ_L)` may not be free

In the inner loop:

```julia
P, N = _clipped_projection_moments( μ[i], μ_L, m_max)
Pm, Nm = _clipped_projection_moments(-μ[i], μ_L, m_max)
```

Under sign flip of μ: `a → -a`, `b → b`, `ψ → π - ψ`. For the half-range
cosine basis used here, sign flip on μ swaps P↔N up to sign on odd m. If
you pin down the exact identity, the second call collapses. Worth ~10
minutes of pencil-and-paper — that's a 2× speedup on the
`_clipped_projection_moments` allocation budget plus correctness
verification (the identity becomes a test invariant).

The `P_m − N_m = a δ_{m0} + (b/2) δ_{m1}` identity already used in the
function points the same direction.

---

## Round 2 — Punch list (recommended priority)

1. ✅ **Issue #1 (Project.toml)** — already fixed.
2. **Rename `*_vsmartmom` helpers** to neutral names. (5 min, big win.)
3. **Loosen `where {FT<:AbstractFloat}` to `<:Real`** — unlock ForwardDiff. (1 line.)
4. **Mutating-variant pattern for `_clipped_projection_moments`** — hoist allocations.
5. **Inline single-use `_canopy_fourier_factor`.**
6. **Strip historical bug comments** before commit; keep convention statements.
7. **Decide on dispatch story** — backend trait vs. three function names. Bigger surgery; do before this is publicly documented.
8. **Make m=0 `compute_Z_matrices` a thin wrapper** around the analytic path; demote brute-force `compute_Z_matrices_aniso(...,m)` to test reference.
9. **Replace `ϖ ≤ 0` throw** with zero-return for the absorbing-leaf limit.
10. **Investigate the `μ → -μ` projection-moment identity** for a 2× speedup and a free correctness invariant.

Items 2, 3, 5, 6, 9 are all sub-10-line edits and worth folding into the
same commit. 4, 7, 8, 10 are follow-ups worth their own commits.

---

## Merging BiLambertian + Specular optical properties

Once the bi-Lambertian fixes land, a natural follow-up is to merge the
two `AbstractCanopyScatteringType` subtypes so a "leaf" can carry both
diffuse and specular components.

### Why the merge is natural

Both produce `(Z⁺⁺, Z⁻⁺)` in the same convention, both take
`(μ, LD, m)`, both already share `AbstractCanopyScatteringType{FT}`.
Real leaves are a *sum* of diffuse + specular contributions
(Verhoef/SAIL, PROSPECT + cuticle Fresnel). Optical properties are
additive at the Z level:

```
Z_leaf = Z_bilambertian + Z_specular
```

…so the cleanest Julian factoring is the same one used for
`CoreScatteringOpticalProperties` in vSmartMOM — overload `+` on the
abstract type:

```julia
struct CompositeCanopyScattering{FT, T<:Tuple} <: AbstractCanopyScatteringType{FT}
    components::T   # (BiLambertianCanopyScattering, SpecularCanopyScattering, ...)
end

Base.:+(a::AbstractCanopyScatteringType, b::AbstractCanopyScatteringType) =
    CompositeCanopyScattering(...)

function compute_Z_matrices(mod::CompositeCanopyScattering, μ, LD, m)
    Zpp, Zmp = compute_Z_matrices(mod.components[1], μ, LD, m)
    for c in Base.tail(mod.components)
        a, b = compute_Z_matrices(c, μ, LD, m)
        Zpp .+= a; Zmp .+= b
    end
    Zpp, Zmp
end
```

Then

```julia
mod = BiLambertianCanopyScattering(R=0.45, T=0.05) +
      SpecularCanopyScattering(nᵣ=1.5, κ=0.2)
```

…and the rest of the canopy code doesn't need to know.

### The constraint

**The closed-form Fourier moments only work for bi-Lambertian.** The
clipped-projection trick relies on the kernel factoring as
`(Ωⁱⁿ·Ωᴸ)(Ωᵒᵘᵗ·Ωᴸ)` so each direction's azimuth integral closes
independently. The specular kernel is localised at the bisector
`Ωᴸ ≈ (Ωⁱⁿ + Ωᵒᵘᵗ) / 2cos(α)` — it doesn't factor that way, and its
Fourier spectrum is broad (sharp peaks → many m). So the specular
component will keep its quadrature-based path indefinitely.

That's *fine* — dispatch handles it: analytic backend for the diffuse
component, brute-force-azimuth backend for the specular component,
summed by the composite. But it's why you don't want to design the API
as "the analytic version" vs "the brute-force version" — it should be
"give me Z up to m_max" and let each component pick its own backend.

### Two cleanups the merge enables

1. **`nQuad` doesn't belong in the scattering structs.** It's an
   integrator parameter, not physics. Both `BiLambertianCanopyScattering`
   and `SpecularCanopyScattering` carry their own `nQuad`, with
   different defaults (30 vs 20). Once components compose, "which
   `nQuad` wins?" becomes a real question. Pull it out as a kwarg on
   `compute_Z_matrices(mod, μ, LD, m; nQuad=...)` or onto a backend
   struct. The bi-Lambertian only needs `nQuad` for the leaf-inclination
   integral now; the specular needs it for outgoing azimuth.
2. **Single-scattering albedo becomes well-defined per component.**
   Bi-Lambertian: `ϖ = R + T`. Specular: `ϖ = ∫ Fresnel·K dΩ` (effective
   specular fraction). Composite: weighted sum. The current `ϖ <= 0`
   throw in the analytic function becomes a non-issue because the
   *composite* ϖ is what matters for normalisation, not any single
   component's.

### Recommended order

1. **Now (this commit)**: land the bi-Lambertian fixes + analytic path.
   Don't touch specular.
2. **Next**: introduce `CompositeCanopyScattering` + `Base.:+`, with
   `compute_Z_matrices(::Composite, ...)` summing components. Pure
   refactor — no behaviour change for callers using a single component.
3. **After**: pull `nQuad` out of the physics structs. This is a
   breaking change for anyone constructing
   `BiLambertianCanopyScattering(R=..., T=..., nQuad=...)` so do it as a
   v0.2 bump.
4. **Later**: tackle the dispatch trait (analytic vs brute-force
   backend) once you actually have a use case for choosing — premature
   otherwise.

The merge is the design pressure that justifies most of the round-2
cleanups (one `compute_Z_matrices` API, no `_vsmartmom` in helper names,
`nQuad` out of the physics struct, `ϖ=0` not an error).
