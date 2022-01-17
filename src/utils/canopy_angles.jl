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

function getSpecularΩ(Ωⁱⁿ::dirVector{FT}, Ωᵒᵘᵗ::dirVector{FT}) where FT
    if Ωⁱⁿ.θ == Ωᵒᵘᵗ.θ && Ωⁱⁿ.ϕ==Ωᵒᵘᵗ.ϕ 
        return Ωⁱⁿ
    end
    θstar  = atan( sqrt(sin(Ωⁱⁿ.θ)^2 + sin(Ωᵒᵘᵗ.θ)^2 - 2sin(Ωⁱⁿ.θ)*sin(Ωᵒᵘᵗ.θ)*cos(Ωⁱⁿ.ϕ - Ωᵒᵘᵗ.ϕ)) / (cos(Ωᵒᵘᵗ.θ) - cos(Ωⁱⁿ.θ)))
    ϕstar  = atan( (sin(Ωⁱⁿ.θ)*sin(Ωⁱⁿ.ϕ) - sin(Ωᵒᵘᵗ.θ)*sin(Ωᵒᵘᵗ.ϕ)) / (sin(Ωⁱⁿ.θ)*cos(Ωⁱⁿ.ϕ) - sin(Ωᵒᵘᵗ.θ)*cos(Ωᵒᵘᵗ.ϕ)))  
    return dirVector(θstar+π/2,ϕstar)
end