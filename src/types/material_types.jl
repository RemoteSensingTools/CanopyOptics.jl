"Abstract Material Type (can be all)"
abstract type AbstractMaterial end

# Subset all water types:
"Abstract water types"
abstract type AbstractWater <: AbstractMaterial end

struct LiquidPureWater <: AbstractWater end

struct LiquidSaltWater <: AbstractWater end

struct PureIce <: AbstractWater end


