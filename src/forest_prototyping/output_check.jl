function check_output_matches_fortran(x::Forest_Scattering_Output)
    
    @assert x.theti == 40

    # hh

    @assert round(x.sighhd1, digits=2) == -16.24
    @assert round(x.sighhd2, digits=2) == -51.46
    @assert round(x.sighhd3, digits=2) == -37.27

    @assert round(x.sighhdr1, digits=3) == -8.708
    @assert round(x.sighhdr2, digits=3) == 8.987
    @assert round(x.sighhdr3, digits=2) == -34.94

    @assert round(x.sighht, digits=3) == 9.073
    @assert round(x.sighhd, digits=2) == -16.24
    @assert round(x.sighhdr, digits=3) == 9.060

    # vv

    @assert round(x.sigvvd1, digits=2) == -10.08
    @assert round(x.sigvvd2, digits=2) == -54.54
    @assert round(x.sigvvd3, digits=2) == -40.02

    @assert round(x.sigvvdr1, digits=3) == -6.647
    @assert round(x.sigvvdr2, digits=3) == 6.000
    @assert round(x.sigvvdr3, digits=2) == -45.73

    @assert round(x.sigvvt, digits=3) == 6.332
    @assert round(x.sigvvd, digits=2) == -10.08
    @assert round(x.sigvvdr, digits=3) == 6.230

    # vh

    @assert round(x.sigvhd1, digits=2) == -17.32
    # @assert round(x.sigvhd2, digits=4) ≈ 0.0
    @assert round(x.sigvhd3, digits=2) == -52.41 

    @assert round(x.sigvhdr1, digits=2) == -12.37
    # @assert round(x.sigvhdr2, digits=4) ≈ 0.0
    @assert round(x.sigvhdr3, digits=2) == -52.01

    @assert round(x.sigvht, digits=2) == -11.16
    @assert round(x.sigvhd, digits=2) == -17.32
    @assert round(x.sigvhdr, digits=2) == -12.37

    # total coherent
    @assert round(x.ghhd, digits=2) == -30.77
    @assert round(x.gvhd, digits=2) == -45.78
    @assert round(x.gvvd, digits=2) == -27.99

    # total incoherent
    @assert round(x.sighhi, digits=3) == 6.075
    @assert round(x.sigvhi, digits=2) == -12.43
    @assert round(x.sigvvi, digits=3) == 3.418

    @assert round(x.sigvhi1, digits=2) == -14.13
    # @assert round(x.sigvhi2, digits=2) ≈ 0
    @assert round(x.sigvhi3, digits=2) == -53.43

    @assert round(x.cofhh, digits=3) == 1.994
    @assert round(x.cofvh, digits=3) == 1.338
    @assert round(x.cofvv, digits=3) == 1.955

    @assert round(x.athc, digits=6) == 0.008110
    @assert round(x.atvc, digits=5) == 0.01734
    @assert round(x.atht, digits=6) == 0.008231
    @assert round(x.atvt, digits=5) == 0.01987

    println("All checks match ✅")

end