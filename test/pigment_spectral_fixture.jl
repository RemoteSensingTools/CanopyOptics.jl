using TOML

function check_pigment_spectra(to_device, synchronize)
    settings=TOML.parsefile(joinpath(@__DIR__, "parameters/pigment_absorption.toml"))
    spectrum=settings["spectrum"]
    native=LeafOpticalProperties()
    for FT in (Float32, Float64)
        optics=LeafOpticalProperties(
            native.λ,
            (
                FT.(getproperty(native, key)) for
                key in (:nᵣ, :Kcab, :Kcar, :Kant, :Kb, :Kw, :Km, :Kp, :Kcbc)
            )...,
        )
        wavelength=FT.(CanopyOptics.ustrip.(native.λ))
        par=findall(x->FT(spectrum["lower_nm"])<=x<FT(spectrum["upper_nm"]), wavelength)
        tolerance=FT(
            settings["validation"][FT==Float32 ? "float32_relative_tolerance" :
                                   "float64_relative_tolerance"],
        )
        photon_energy=FT(spectrum["planck_J_s"]) *
                      FT(spectrum["light_speed_m_s"]) *
                      FT(spectrum["avogadro_mol_inverse"])/FT(
            spectrum["micromoles_per_mole"],
        ) ./ (wavelength[par] .* FT(spectrum["nm_to_m"]))
        for raw in settings["leaf"]
            names=(:N, :Ccab, :Ccar, :Canth, :Cbrown, :Cw, :Cm, :Cprot, :Ccbc)
            leaf=LeafProspectProProperties{FT}(;
                (key=>FT(raw[string(key)]) for key in names)...,
            )
            @test iszero(leaf.Cm)
            budget=prospect_budget(leaf, optics)
            products=@inferred prospect_constituent_coefficients(leaf, optics)
            @test all(x->eltype(x)===FT, values(products))
            @test products.chlorophyll==leaf.Ccab*optics.Kcab
            @test products.carotenoids==leaf.Ccar*optics.Kcar
            @test products.protein==leaf.Cprot*optics.Kp
            @test products.carbon_constituents==leaf.Ccbc*optics.Kcbc
            per_bin=[
                NamedTuple{keys(products)}(map(x->x[i], values(products))) for i in par
            ]
            absorbed=budget.absorptance[par] .* FT(spectrum["incident_photons_per_bin"])
            result=partition_absorption.(
                Ref(HomogeneousAbsorption()),
                to_device(absorbed),
                to_device(per_bin),
            )
            synchronize()
            collected=Array(result)
            @test all(absorption_valid, collected)
            @test all(x->x.status===AbsorptionPartitionValid, collected)
            # Check both photon closure and its wavelength-resolved energy weighting.
            independent_photon_sum=sum(BigFloat(q) for q in absorbed)
            independent_energy_sum=sum(
                BigFloat(q)*BigFloat(e) for (q, e) in zip(absorbed, photon_energy)
            )
            species_photons=sum(
                BigFloat(q) for row in collected for q in values(row.constituents)
            )
            species_energy=sum(
                BigFloat(q)*BigFloat(e) for (row, e) in zip(collected, photon_energy) for
                q in values(row.constituents)
            )
            @test isapprox(species_photons, independent_photon_sum; rtol = tolerance)
            @test isapprox(species_energy, independent_energy_sum; rtol = tolerance)
            @test maximum(
                abs.(budget.transmission .+ budget.reflectance .+ budget.absorptance .- 1),
            )<=12eps(FT)
            previous=zero(FT)
            for eta in FT.(settings["selection"]["carotenoid_transfer_fractions"])
                selected=antenna_excitation.(Ref(CarotenoidAntenna(eta)), result)
                psii=photosystem_excitation.(
                    selected,
                    FT(settings["selection"]["photosystem_fraction"]),
                )
                synchronize()
                total=sum(selected)
                @test zero(FT)<=sum(psii)<=total<=sum(to_device(absorbed))
                @test previous<=total
                previous=total
                for (a, row) in zip(Array(selected), collected)
                    reference=BigFloat(row.constituents.chlorophyll)+BigFloat(eta)*BigFloat(
                        row.constituents.carotenoids,
                    )
                    @test isapprox(a, FT(reference); rtol = tolerance)
                end
            end
            chl=antenna_excitation.(Ref(ChlorophyllOnly()), result)
            eta0=antenna_excitation.(Ref(CarotenoidAntenna(zero(FT))), result)
            @test Array(chl)==Array(eta0)
        end
    end
end
