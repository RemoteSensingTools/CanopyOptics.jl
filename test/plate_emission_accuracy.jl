# Regenerate the scalar accuracy record after running runtests.jl and CUDA tests.
using CanopyOptics, TOML, SHA, Dates
include("plate_emission_reference.jl")
include("plate_emission_batches.jl")
root=dirname(@__DIR__)
records=Dict{String,Any}()
for FT in (Float32, Float64)
    escapes, births, expected=plate_emission_batch_cases(FT)
    # Include the direct/Taylor transition and the excitation moment switch.
    threshold=sqrt(sqrt(eps(FT)))
    for a in
        FT.((prevfloat(FT(1)), 1, nextfloat(FT(1)), prevfloat(FT(8)), 8, nextfloat(FT(8)))),
        b in (prevfloat(threshold), threshold, nextfloat(threshold)),
        t in FT.((1, 0.6, 1e-5))

        push!(escapes, RegularPlateEscape(b, t))
        push!(births, ExponentialBirthProfile(a, one(FT), FT(0.3)))
        push!(expected, map(FT, plate_emission_reference(a, b, t, one(FT), FT(0.3))))
    end
    relative=Dict(string(f)=>0.0 for f in (:front, :back, :reabsorbed))
    absolute=copy(relative)
    worst=Dict{String,Any}()
    closure=0.0
    for (escape, birth, ref) in zip(escapes, births, expected)
        value=leaf_emission_fates(escape, birth)
        closure=max(closure, Float64(abs(sum(values(value))-one(FT))))
        for field in (:front, :back, :reabsorbed)
            v, r=getproperty(value, field), getproperty(ref, field)
            error=abs(BigFloat(v)-BigFloat(r))
            absolute[string(field)]=max(absolute[string(field)], Float64(error))
            if abs(r)>=floatmin(FT)
                err=Float64(error/abs(BigFloat(r)))
                if err>relative[string(field)]
                    relative[string(field)]=err
                    worst[string(field)]=Dict(
                        "excitation_depth"=>Float64(birth.excitation_depth),
                        "emission_depth"=>Float64(escape.emission_depth),
                        "interface_transmission"=>Float64(escape.interface_transmission),
                    )
                end
            end
        end
    end
    records[string(FT)]=Dict(
        "cases"=>length(escapes),
        "max_relative_error_normal_outputs"=>relative,
        "max_absolute_error"=>absolute,
        "max_absolute_photon_closure"=>closure,
        "worst_relative_case"=>worst,
    )
end
files=[
    "src/CanopyOptics.jl",
    "src/leaf_optics/plate_emission.jl",
    "test/test_plate_emission.jl",
    "test/plate_emission_reference.jl",
    "test/plate_emission_batches.jl",
    "test/plate_emission_gpu.jl",
    "test/plate_emission_accuracy.jl",
    "test/runtests.jl",
    "test/reference/generate_plate_emission_reference.py",
    "test/reference/plate_emission_reference.csv",
]
report=Dict(
    "generated_utc"=>string(now(UTC)),
    "julia_version"=>string(VERSION),
    "reference"=>"2048-bit independent analytic boundary balances; separate 90-digit quadrature/path fixtures",
    "accuracy"=>records,
    "sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(root, f)))) for f in files),
)
logdir=only(ARGS)
report["validation_logs"]=Dict(
    name=>Dict(
        "sha256"=>bytes2hex(sha256(read(joinpath(logdir, name)))),
        "text"=>read(joinpath(logdir, name), String),
    ) for name in ("full.log", "gpu.log", "docs.log")
)
output=joinpath(root, "docs", "src", "data", "plate_emission_validation.toml")
mkpath(dirname(output))
open(output, "w") do io
    TOML.print(io, report; sorted = true)
end
println(output)
TOML.print(stdout, records; sorted = true)
