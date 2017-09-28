using FinEtools
using FinEtools.AlgoDeforLinearModule

function calculateerrors!(coarse, fine, finest)
    # Load coarse-mesh data
    fensc = coarse["fens"]
    fesc = coarse["fes"]
    stressfieldsc = coarse["stressfields"]
    uc = coarse["u"]

    # Load fine-mesh data
    fensf = fine["fens"]
    fesf = fine["fes"]
    femmf = fine["femm"]
    stressfieldsf = fine["stressfields"]
    geomf = fine["geom"]
    uf = fine["u"]

    # Load finest-mesh data
    fensfst = finest["fens"]
    fesfst = finest["fes"]
    femmfst = finest["femm"]
    stressfieldsfst = finest["stressfields"]
    geomfst = finest["geom"]
    ufst = finest["u"]
    tolerance = fine["tolerance"]

    # (a) Displacement
    # Transfer the result from the coarse mesh to the fine-mesh
    fieldfst = ufst
    fieldf = uf
    fieldc = uc
    fieldt = NodalField(zeros(size(fieldf.values)))

    fieldt = transferfield!(fieldt, fensf, fesf, fieldc, fensc, fesc, tolerance)

    diffff = NodalField(fieldf.values - fieldt.values)
    geom = NodalField(fensf.xyz)
    errornorm = integratefieldfunction(femmf, geomf, diffff, (x, v) -> norm(v)^2, 0.0)
    errornorm = sqrt(errornorm)
    solnorm = integratefieldfunction(femmfst, geomfst, fieldfst, (x, v) -> norm(v)^2, 0.0)
    solnorm = sqrt(solnorm)
    println("|error|/|ref| = $(errornorm/solnorm)")
    coarse["displacement_errornorm"] = errornorm
    coarse["displacement_solnorm"] = solnorm


    # (b) Stresses
    # Transfer the result from the coarse mesh to the fine-mesh
    stress = FDataDict[]
    i = 1 # If there should be more fields  in stressfieldsf, this would have to change
    fieldfst = stressfieldsfst[i]
    fieldf = stressfieldsf[i]
    fieldc = stressfieldsc[i]
    fieldt = deepcopy(fieldf)

    fieldt = transferfield!(fieldt, fensf, fesf, fieldc, fensc, fesc, tolerance)

    diffff = deepcopy(fieldf)
    diffff.values[:] = fieldf.values - fieldt.values
    geom = NodalField(fensf.xyz)
    errornorm = integratefieldfunction(femmf, geomf, diffff, (x, v) -> norm(v)^2, 0.0)
    errornorm = sqrt(errornorm)
    solnorm = integratefieldfunction(femmfst, geomfst, fieldfst, (x, v) -> norm(v)^2, 0.0)
    solnorm = sqrt(solnorm)
    println("|error|/|ref| = $(errornorm/solnorm)")
    push!(stress, FDataDict())
    coarse["stress_errornorm"] = errornorm
    coarse["stress_solnorm"] = solnorm

    return coarse
end


File = "fiber_reinf_cant_iso_stresses_MST10"
convergencestudy = open(File * ".jls", "r") do file
    deserialize(file)
end

elementsizes = [s["elementsize"] for s in convergencestudy]
println("elementsizes = $(elementsizes)")

displacementerrornorm = FFlt[]
displacementsolnorm = FFlt[]
stresserrornorm = FFlt[]
stresssolnorm = FFlt[]
for i = 1:(length(convergencestudy) - 1)
    calculateerrors!(convergencestudy[i], convergencestudy[i + 1], convergencestudy[end])
    push!(displacementerrornorm, convergencestudy[i]["displacement_errornorm"])
    push!(displacementsolnorm, convergencestudy[i]["displacement_solnorm"])
    push!(stresserrornorm, convergencestudy[i]["stress_errornorm"])
    push!(stresssolnorm, convergencestudy[i]["stress_solnorm"])
end

open(File * ".jls", "w") do file
    serialize(file, convergencestudy)
end

println("displacementerrornorm = $(displacementerrornorm)")
println("displacementsolnorm = $(displacementsolnorm)")
println("stresserrornorm = $(stresserrornorm)")
println("stresssolnorm = $(stresssolnorm)")

csvFile = File * "_errors" * ".CSV"
savecsv(File * "_errors" * ".CSV",
    elementsizes=vec(elementsizes[1:end-1]),
    elementsizes2=vec(elementsizes[1:end-1].^2),
    u_normerror=vec(displacementerrornorm ./ displacementsolnorm),
    s_normerror=vec(stresserrornorm ./ stresssolnorm),
    )

@async run(`"paraview.exe" $csvFile`)
