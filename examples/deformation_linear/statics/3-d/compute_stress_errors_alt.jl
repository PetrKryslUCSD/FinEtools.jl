using FinEtools
using FinEtools.AlgoDeforLinearModule
using JLD

function calculateerrors(coarse, benchmark)
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


    # Load benchmark-mesh data
    File = benchmark
    fensb = jldopen(File, "r") do file
        read(file, "fens")
    end
    fesb = jldopen(File, "r") do file
        read(file, "fes")
    end
    stressfieldsb = jldopen(File, "r") do file
        read(file, "stressfields")
    end
    geomb = jldopen(File, "r") do file
        read(file, "geom")
    end
    ub = jldopen(File, "r") do file
        read(file, "u")
    end
    integrationrule = jldopen(File, "r") do file
        read(file, "integrationrule")
    end
    tolerance = jldopen(File, "r") do file
        read(file, "tolerance")
    end

    # (a) Displacement
    # Transfer the result from the coarse mesh to the fine-mesh
    displacement = FDataDict()
    fieldb = ub
    fieldc = uc
    fieldt = NodalField(zeros(size(fieldb.values)))

    fieldt = transferfield!(fieldt, fensb, fesb, fieldc, fensc, fesc, tolerance)

    diffff = NodalField(fieldb.values - fieldt.values)
    geomb = NodalField(fensb.xyz)
    femmb = FEMMBase(GeoD(fesb, integrationrule))
    normerror = integratefieldfunction(femmb, geomb, diffff, (x, v) -> norm(v)^2, 0.0)
    normerror = sqrt(normerror)
    normref = integratefieldfunction(femmb, geomb, fieldb, (x, v) -> norm(v)^2, 0.0)
    normref = sqrt(normref)
    println("|error|/|ref| = $(normerror/normref)")
    displacement["normerror"] = normerror
    displacement["normref"] = normref

    # (b) Stresses
    # Transfer the result from the coarse mesh to the fine-mesh
    stress = FDataDict[]
    for i = 1:length(stressfieldsb)
        fieldb = stressfieldsb[i]
        fieldc = stressfieldsc[i]
        fieldt = deepcopy(fieldb)

        fieldt = transferfield!(fieldt, fensb, fesb, fieldc, fensc, fesc, tolerance)

        diffff = deepcopy(fieldb)
        diffff.values[:] = fieldb.values - fieldt.values
        geomb = NodalField(fensb.xyz)
        femmb = FEMMBase(GeoD(fesb, integrationrule))
        normerror = integratefieldfunction(femmb, geomb, diffff, (x, v) -> norm(v)^2, 0.0)
        normerror = sqrt(normerror)
        normref = integratefieldfunction(femmb, geomb, fieldb, (x, v) -> norm(v)^2, 0.0)
        normref = sqrt(normref)
        println("|error|/|ref| = $(normerror/normref)")
        push!(stress, FDataDict())
        stress[i]["normerror"] = normerror
        stress[i]["normref"] = normref
    end

    return displacement, stress

end

# displacementerror = []
# sxxerror = []
# sxzerror = []
# ns =  [2 4 8 16]
# for n = ns
#     coarse = "fiber_reinf_cant_iso_stressesn=$(n).jld"
#     fine = "fiber_reinf_cant_iso_stressesn=$(2*n).jld"
#     displacement, stress = calculateerrors(coarse, fine)
#     push!(displacementerror, displacement["normerror"]/displacement["normref"])
#     push!(sxxerror, stress[1]["normerror"]/stress[1]["normref"])
#     push!(sxzerror, stress[5]["normerror"]/stress[5]["normref"])
# end
#
# File = "fiber_reinf_cant_iso_stresses_errors" * ".CSV"
# savecsv(File, ns=vec(ns),
#     displacementerror=vec(displacementerror),
#     sxxerror=vec(sxxerror),
#     sxzerror=vec(sxzerror),
#     ns1=vec(1.0 ./ ns),
#     ns2=vec(1.0 ./ ns.^2))

displacementerror = []
stresserror = []
ns =  [1 2 4 8]
for n = ns
    coarse = "fiber_reinf_cant_iso_stressesn=$(n).jld"
    benchmark = "fiber_reinf_cant_iso_stressesn=16.jld"
    displacement, stress = calculateerrors(coarse, benchmark)
    push!(displacementerror, displacement["normerror"]/displacement["normref"])
    push!(stresserror, stress[1]["normerror"]/stress[1]["normref"])
end
#
# println("stress[1]["normref"] = $(stress[1]["normref"])")


File = "fiber_reinf_cant_iso_stresses_errors" * ".CSV"
savecsv(File, ns=vec(ns),
    displacementerror=vec(displacementerror),
    stresserror=vec(stresserror),
    ns1=vec(1.0 ./ ns),
    ns2=vec(1.0 ./ ns.^2))

@async run(`"paraview.exe" $File`)
