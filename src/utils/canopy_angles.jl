"Get Beta Distribution parameters using standard tabulated values θₗ,s, Eq. 2.9 and 2.10 in Bonan"
function βparameters(θₗ,s)
    α = (1 - (s^2+θₗ^2) / (θₗ*π/2) ) / ( ((s^2+θₗ^2) / θₗ^2) - 1 )
    β = ( (π/2)/θₗ -1 ) * α
    return α,β
end

"Eq 46 in Shultis and Myneni"
function H(μ::FT,μₗ::FT) where FT
    # cot(θ) * cot(θₗ)
    r = cot(acos(μ)) * cot(acos(μₗ));
    if r > 1
        return μ * μₗ
    elseif r < -1
        return FT(0)
    else
        # Eq 47
        ϕ = acos(-r);
        return FT(1/π) * (μ*μₗ*ϕ + sqrt(FT(1)-μ^2) * sqrt(1-μₗ^2) * sin(ϕ)) 
    end
end

"Eq 45 in Shultis and Myneni, fixed grid in μ for both μ and μ' "
function compute_Ψ(μ::Array{FT,1}, μₗ::FT) where FT
    H⁺ = H.(μ,μₗ)
    H⁻ = H.(-μ,μₗ)
    Ψ⁺ = H⁺ * H⁺' + H⁻ * H⁻'
    Ψ⁻ = H⁺ * H⁻' + H⁻ * H⁺'
    return Ψ⁺, Ψ⁻
end

"Eq 45 in Shultis and Myneni, fixed grid in μ for both μ and μ' "
function compute_Ψ(μ::Array{FT,1},μꜛ::Array{FT,1}, μₗ::FT) where FT
    Ψ⁺ = H.(μ,μₗ) * H.(μꜛ,μₗ)'  + H.(-μ,μₗ) * H.(-μꜛ,μₗ)'
    Ψ⁻ = H.(μ,μₗ) * H.(-μꜛ,μₗ)' + H.(-μ,μₗ) * H.(+μꜛ,μₗ)'
    return Ψ⁺, Ψ⁻
end

" Eq 37 in Shultis and Myneni "
function scattering_model_lambertian(Ωⁱⁿ::dirVector, Ωᵒᵘᵗ::dirVector, Ωᴸ::dirVector, r, t)
    α = (Ωⁱⁿ ⋅ Ωᴸ) * (Ωᵒᵘᵗ ⋅ Ωᴸ)
    if a<0
        return r * abs(Ωᵒᵘᵗ ⋅ Ωᴸ) / π
    else
        return t * abs(Ωᵒᵘᵗ ⋅ Ωᴸ) / π
    end
end

#See Eq 18, Canopy RT book
function getSpecularΩ(Ωⁱⁿ::dirVector{FT}, Ωᵒᵘᵗ::dirVector{FT}) where FT
    sa = Ωⁱⁿ ⋅ Ωᵒᵘᵗ 
    sa > 1 ? sa = FT(1) : nothing
    
    α = acos(sa)/2
    # Relative azimuth angle:
    ϕ = (Ωᵒᵘᵗ.ϕ - Ωⁱⁿ.ϕ) 
    # Leaf polar angle:
    θₗ = acos((cos(Ωⁱⁿ.θ) + cos(Ωᵒᵘᵗ.θ))/2cos(α));
    
    t = (sin(Ωⁱⁿ.θ) + sin(Ωᵒᵘᵗ.θ) * cos(ϕ))/(2cos(α)*sin(θₗ))
    if isnan(t)
        t = FT(1.0)
    end
    # Leaf azimuth angle:
    ϕₗ = acos(max(FT(-1),min(t,FT(1))))
    ϕ > π   ? ϕₗ = 2π-ϕₗ-Ωⁱⁿ.ϕ : ϕₗ = ϕₗ+Ωⁱⁿ.ϕ
    #θₗ > π/2 ? θₗ = θₗ - π/2  : nothing 
    out = dirVector(θₗ,ϕₗ)
    #ϕ > π   ? out = dirVector(θₗ,2π-ϕₗ-Ωⁱⁿ.ϕ) : out = dirVector(θₗ,ϕₗ+Ωⁱⁿ.ϕ)
    #@show (Ωⁱⁿ ⋅ out) , (Ωᵒᵘᵗ ⋅ out)
    #@show abs(((Ωⁱⁿ ⋅ out) - (Ωᵒᵘᵗ ⋅ out)) / (Ωᵒᵘᵗ ⋅ out)) * 100
    #@show θₗ , ϕₗ 
    @assert abs(((Ωⁱⁿ ⋅ out) - (Ωᵒᵘᵗ ⋅ out)) / (Ωᵒᵘᵗ ⋅ out)) * 100 < 1 "Leaf angle wrong, $Ωⁱⁿ, $Ωᵒᵘᵗ, $out"
    return out 

end

function getSpecularΩ(Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT}) where FT
    sa = Ωⁱⁿ ⋅ Ωᵒᵘᵗ 
    #@show sa
    sa > 1 ? sa = FT(1) : nothing
    # Scattering angle
    α = acos(sa)/FT(2)
    
    # Relative azimuth angle:
    ϕ = (Ωᵒᵘᵗ.ϕ - Ωⁱⁿ.ϕ) 
    # Leaf polar angle:
    μₗ = ((Ωⁱⁿ.μ + Ωᵒᵘᵗ.μ)/2cos(α));
    
    t = (sqrt(1-Ωⁱⁿ.μ^2) + sqrt(1-Ωᵒᵘᵗ.μ^2) * cos(ϕ))/(2cos(α)*sqrt(1-μₗ^2))
    #@show t, sqrt(1-Ωⁱⁿ.μ^2) + sqrt(1-Ωᵒᵘᵗ.μ)
    if isnan(t)
        t = FT(1.0)
    end
    # Leaf azimuth angle:
    ϕₗ = acos(max(FT(-1),min(t,FT(1))))
    ϕ > π   ? ϕₗ = 2π-ϕₗ-Ωⁱⁿ.ϕ : ϕₗ = ϕₗ+Ωⁱⁿ.ϕ
    out = dirVector_μ(μₗ,ϕₗ)
    @assert Ωⁱⁿ ⋅ out ≈  Ωᵒᵘᵗ ⋅ out "Leaf angle wrong in specular reflection, $Ωⁱⁿ, $Ωᵒᵘᵗ, $out"
    return out, α
end