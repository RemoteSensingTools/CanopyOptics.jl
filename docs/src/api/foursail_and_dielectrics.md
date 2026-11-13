# FourSAIL & Microwave Dielectrics

This page documents the `FourSAIL` four-stream turbid-medium canopy radiative-transfer
solver, complex microwave relative permittivity models ``\varepsilon(T, f)`` ([`dielectric`](@ref)), and
deprecated compatibility aliases.

```@contents
Pages = ["foursail_and_dielectrics.md"]
Depth = 2
```

## FourSAIL Canopy Radiative Transfer

```@docs
FourSAILGeometry
FourSAILGeometrySet
FourSAILResult
foursail
foursail!
```

## Microwave Dielectric Permittivity (`ε(T, f)`)

```@docs
AbstractMaterial
AbstractWater
AbstractSoil
AbstractVegetation
LiquidPureWater
LiquidSaltWater
PureIce
SoilMW
LeafUlabyElRayes1987
dielectric
```

## Deprecated Compatibility Aliases

```@docs
createLeafOpticalStruct
createMultiPigmentOpticalStruct
dirVector
dirVector_μ
```
