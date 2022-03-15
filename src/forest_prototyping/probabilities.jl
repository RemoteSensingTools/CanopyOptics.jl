# θ:      input value
# ntype:  distribution type
# param:  distribution parameter
function prob(θ, ntype::Integer, param=0.0)
    @assert 0 <= θ <= π "θ must be between 0 and π"
    θ_d = θ*180.0/pi
    if ntype == 1
        return (0 < θ < π/2 ? 0.5*sin(θ) : 0.0)
    elseif ntype == 2
        return (0 < θ < π/2 ? 8/(3π) * sin(θ)^4 : 0.0)
    # Unfinished (3)
    elseif ntype == 3
        global nth
        th0r = param*pi/180.
        dth = π/nth
        return (θ <= th0r < θ + dth ? 1.0/dth : 0.0)
    elseif ntype == 4
        dthplr = param * π/180
        return (0 < θ < dthplr ? 1.0/dthplr : 0.0)
    elseif ntype == 5
        dthpdr = param * π/180
        return ((π/2 - dthpdr) < θ < (π/2 + dthpdr) ? 1.0/(2.0*dthpdr) : 0.0)
    elseif ntype == 6
        return (θ < π/2 ? (16.0/3.0*π)*(sin(2.0*θ))^4 : 0.0)
    elseif ntype == 7
        return (π/4 < θ < π/2 ? abs(2 * sin(4.0 * θ)) : 0.0) 
    elseif ntype == 8
        return (16.0/(5.0*pi))*(cos(θ))^6
    elseif ntype == 9
        return param * exp(-param * θ)
    elseif ntype == 10
        thd = θ*180.0/pi
        if (0 <=thd <= 10.0)
            return 0.015
        elseif (10.0 < thd <= 20.0)
            return 0.020 
        elseif (20.0 < thd <= 30.0)
            return 0.015
        elseif (30.0 < thd <= 40.0)
            return 0.03
        elseif (40.0 < thd <= 50.0)
            return 0.14
        elseif (50.0 < thd <= 60.0)
            return 0.25
        elseif (60.0 < thd <= 70.0)
            return 0.22
        elseif (70.0 < thd <= 80.0)
            return 0.14
        elseif (80.0 < thd <= 90.0)
            return 0.17
        end
    elseif ntype == 11
        return 2 * (1 + cos(2*θ))/pi
    elseif ntype == 12
        return 2 * (1 - cos(2*θ))/pi
    elseif ntype == 13
        return 2 * (1 - cos(4*θ))/pi
    elseif ntype == 14
        return 2 * (1 + cos(4*θ))/pi
    elseif ntype == 15
        return 2/π
    elseif ntype == 16
        return sin(θ)
    elseif ntype == 17
        if (0 < θ_d <= 40.0)
            return 0.0
        elseif (40.0 < θ_d <= 60.0)
            return 0.053
        elseif (60.0 < θ_d <= 80.0)
            return 0.149
        elseif (80.0 < θ_d <= 100.0)
            return 0.255
        elseif (100.0 < θ_d <= 120.0)
            return 0.416
        elseif (120.0 < θ_d <= 140.0)
            return 0.114
        elseif (140.0 < θ_d)
            return 0.0
        end
    end
end