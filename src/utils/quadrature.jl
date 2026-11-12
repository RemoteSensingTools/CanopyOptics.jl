"""
    gauleg(n, xmin, xmax; norm = false)

Compute `n`-point Gauss–Legendre quadrature nodes and weights mapped onto
`[xmin, xmax]`. When `norm = true`, weights are normalized to sum to `1`.

# Examples
```jldoctest
julia> x, w = CanopyOptics.gauleg(4, 0.0, 1.0);

julia> isapprox(sum(w), 1.0; atol = 1e-15)
true
```
"""
function gauleg(n::Integer, xmin_in::Real, xmax_in::Real; norm::Bool = false)
    FT = float(promote_type(typeof(xmin_in), typeof(xmax_in)))
    xmin = FT(xmin_in)
    xmax = FT(xmax_in)
    ξ_raw, w_raw = gausslegendre(Int(n))
    half_span = (xmax - xmin) / FT(2)
    mid = (xmin + xmax) / FT(2)
    ξ = FT.(half_span .* ξ_raw .+ mid)
    w = norm ? FT.(w_raw ./ sum(w_raw)) : FT.(w_raw .* half_span)
    return ξ, w
end
