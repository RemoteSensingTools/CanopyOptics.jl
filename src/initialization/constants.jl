# Empirical polynomial and relaxation coefficients for Ulaby & Long (2014)
# saline and pure liquid water permittivity models (Chapter 4, Eqs. 4.19–4.21).
const _TEMP_POLY_COEFFS   = (2.903602, 8.607e-2, 4.738817e-4, -2.991e-6, 4.3041e-9)
const _SAL_POLY1_COEFFS   = (37.5109, 5.45216, 0.014409)
const _SAL_POLY2_COEFFS   = (1004.75, 182.283, 1.0)
const _ALPHA0_NOM_COEFFS  = (6.9431, 3.2841, -0.099486)
const _ALPHA0_DEN_COEFFS  = (84.85, 69.024, 1.0)
const _ALPHA1_COEFFS      = (49.843, -0.2276, 0.00198)

const _ULABY_WATER_AU = (
    0.46606917e-2, -0.26087876e-4, -0.63926782e-5,
    0.63000075e1,   0.26242021e-2, -0.42984155e-2,
    0.34414691e-4,  0.17667420e-3, -0.20491560e-6,
    0.58366888e3,   0.12634992e3,   0.69227972e-4,
    0.38957681e-6,  0.30742330e3,   0.12634992e3,
    0.37245044e1,   0.92609781e-2, -0.26093754e-1,
)
