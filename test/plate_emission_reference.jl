# Independent high-precision integration of the exponential birth density and
# solution of the two boundary balance equations. No production helper calls.
function plate_emission_reference(a, b, t, f, g)
    setprecision(BigFloat, 2048) do
        a, b, t, f, g = BigFloat.((a, b, t, f, g))
        L(x) = iszero(x) ? one(x) : (1-exp(-x))/x
        S = L(a+b)
        D = a == b ? exp(-a) : (exp(-b)-exp(-a))/(a-b)
        normalization = (f+g)*L(a)
        first_front = (f*S+g*D)/(2normalization)
        first_back = (g*S+f*D)/(2normalization)
        r = 1-t
        q = r*exp(-b)
        # The inverse of [1 -q; -q 1], evaluated at 2048 bits.
        front = t*(first_front+q*first_back)/(1-q*q)
        back = t*(first_back+q*first_front)/(1-q*q)
        return (; front, back, reabsorbed = 1-front-back)
    end
end
