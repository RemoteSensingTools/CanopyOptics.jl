function check_output_matches_fortran(x::Forest_Scattering_Output)
    
    @test x.theti == 40

    # hh

    @test round(x.sighhd1, digits=2) == -16.24
    @test round(x.sighhd2, digits=2) == -51.46
    @test round(x.sighhd3, digits=2) == -37.27

    @test round(x.sighhdr1, digits=3) == -8.708
    @test round(x.sighhdr2, digits=3) == 8.987
    @test round(x.sighhdr3, digits=2) == -34.94

    @test round(x.sighht, digits=3) == 9.073
    @test round(x.sighhd, digits=2) == -16.24
    @test round(x.sighhdr, digits=3) == 9.060

    # vv

    @test round(x.sigvvd1, digits=2) == -10.08
    @test round(x.sigvvd2, digits=2) == -54.54
    @test round(x.sigvvd3, digits=2) == -40.02

    @test round(x.sigvvdr1, digits=3) == -6.647
    @test round(x.sigvvdr2, digits=3) == 6.000
    @test round(x.sigvvdr3, digits=2) == -45.73

    @test round(x.sigvvt, digits=3) == 6.332
    @test round(x.sigvvd, digits=2) == -10.08
    @test round(x.sigvvdr, digits=3) == 6.230

    # vh

    @test round(x.sigvhd1, digits=2) == -17.32
    # @test round(x.sigvhd2, digits=4) ≈ 0.0
    @test round(x.sigvhd3, digits=2) == -52.41 

    @test round(x.sigvhdr1, digits=2) == -12.37
    # @test round(x.sigvhdr2, digits=4) ≈ 0.0
    @test round(x.sigvhdr3, digits=2) == -52.01

    @test round(x.sigvht, digits=2) == -11.16
    @test round(x.sigvhd, digits=2) == -17.32
    @test round(x.sigvhdr, digits=2) == -12.37

    # total coherent
    @test round(x.ghhd, digits=2) == -30.77
    @test round(x.gvhd, digits=2) == -45.78
    @test round(x.gvvd, digits=2) == -27.99

    # total incoherent
    @test round(x.sighhi, digits=3) == 6.075
    @test round(x.sigvhi, digits=2) == -12.43
    @test round(x.sigvvi, digits=3) == 3.418

    @test round(x.sigvhi1, digits=2) == -14.13
    # @test round(x.sigvhi2, digits=2) ≈ 0
    @test round(x.sigvhi3, digits=2) == -53.43

    @test round(x.cofhh, digits=3) == 1.994
    @test round(x.cofvh, digits=3) == 1.338
    @test round(x.cofvv, digits=3) == 1.955

    @test round(x.athc, digits=6) == 0.008110
    @test round(x.atvc, digits=5) == 0.01734
    @test round(x.atht, digits=6) == 0.008231
    @test round(x.atvt, digits=5) == 0.01987

    println("All checks match ✅")

end