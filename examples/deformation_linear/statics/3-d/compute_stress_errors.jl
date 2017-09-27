using FinEtools
using FinEtools.AlgoDeforLinearModule
using JLD

function calculateerrors(coarse, fine)
    # Load coarse-mesh data
    File = coarse
    fensc = jldopen(File, "r") do file
        read(file, "fens")
    end
    fesc = jldopen(File, "r") do file
        read(file, "fes")
    end
    stressfieldsc = jldopen(File, "r") do file
        read(file, "stressfields")
    end
    uc = jldopen(File, "r") do file
        read(file, "u")
    end

    # Load fine-mesh data
    File = fine
    fensf = jldopen(File, "r") do file
        read(file, "fens")
    end
    fesf = jldopen(File, "r") do file
        read(file, "fes")
    end
    stressfieldsf = jldopen(File, "r") do file
        read(file, "stressfields")
    end
    geomf = jldopen(File, "r") do file
        read(file, "geom")
    end
    uf = jldopen(File, "r") do file
        read(file, "u")
    end
    # femmf = jldopen(File, "r") do file
    #     read(file, "femm")
    # end
    integrationrule = jldopen(File, "r") do file
        read(file, "integrationrule")
    end
    tolerance = jldopen(File, "r") do file
        read(file, "tolerance")
    end

    # (a) Displacement
    # Transfer the result from the coarse mesh to the fine-mesh
    displacement = FDataDict()
    fieldf = uf
    fieldc = uc
    fieldt = NodalField(zeros(size(fieldf.values)))

    fieldt = transferfield!(fieldt, fensf, fesf, ieldc, fensc, fesc, tolerance)

    diffff = NodalField(fieldf.values - fieldt.values)
    geom = NodalField(fensf.xyz)
    femmf = FEMMBase(GeoD(fesf, integrationrule))
    normerror = integratefieldfunction(femmf, geomf, diffff, (x, v) -> norm(v), 0.0)
    normref = integratefieldfunction(femmf, geomf, fieldf, (x, v) -> norm(v), 0.0)
    println("|error|/|ref| = $(normerror/normref)")
    displacement["normerror"] = normerror
    displacement["normref"] = normref

    # (b) Stresses
    # Transfer the result from the coarse mesh to the fine-mesh
    stress = FDataDict[]
    for i = 1:length(stressfieldsf)
        fieldf = stressfieldsf[i]
        fieldc = stressfieldsc[i]
        fieldt = NodalField(zeros(size(fieldf.values)))

        fieldt = transferfield!(fieldt, fensf, fesf, fieldc, fensc, fesc, tolerance)

        diffff = NodalField(fieldf.values - fieldt.values)
        geom = NodalField(fensf.xyz)
        femmf = FEMMBase(GeoD(fesf, integrationrule))
        normerror = integratefieldfunction(femmf, geomf, diffff, (x, v) -> norm(v), 0.0)
        normref = integratefieldfunction(femmf, geomf, fieldf, (x, v) -> norm(v), 0.0)
        println("|error|/|ref| = $(normerror/normref)")
        push!(stress, FDataDict())
        stress[i]["normerror"] = normerror
        stress[i]["normref"] = normref
    end

    return displacement, stress

end

displacementerror = []
sxxerror = []
sxzerror = []
ns =  [2 4 8 16]
for n = ns
    coarse = "fiber_reinf_cant_iso_stressesn=$(n).jld"
    fine = "fiber_reinf_cant_iso_stressesn=$(2*n).jld"
    displacement, stress = calculateerrors(coarse, fine)
    push!(displacementerror, displacement["normerror"]/displacement["normref"])
    push!(sxxerror, stress[1]["normerror"]/stress[1]["normref"])
    push!(sxzerror, stress[5]["normerror"]/stress[5]["normref"])
end

File = "fiber_reinf_cant_iso_stresses_errors" * ".CSV"
savecsv(File, ns=vec(ns),
    displacementerror=vec(displacementerror),
    sxxerror=vec(sxxerror),
    sxzerror=vec(sxzerror),
    ns1=vec(1.0 ./ ns),
    ns2=vec(1.0 ./ ns.^2))

@async run(`"paraview.exe" $File`)
