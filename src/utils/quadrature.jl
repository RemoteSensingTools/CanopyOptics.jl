function gauleg(n::Int, xmin::FT, xmax::FT; norm=false) where FT
    ξ, w = gausslegendre(n)
    ξ = (xmax - xmin) / FT(2) * ξ .+ (xmin + xmax) / FT(2)
    norm ? w /= sum(w) : w *= (xmax - xmin) / FT(2)
    return FT.(ξ), FT.(w)
end