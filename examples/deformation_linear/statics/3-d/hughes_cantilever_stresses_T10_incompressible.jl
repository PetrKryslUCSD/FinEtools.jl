using FinEtools
using FinEtools.AlgoDeforLinearModule
using ComputeErrorsModule

elementtag = "T10"
println("""
Cantilever example.  Hughes 1987. Element: $(elementtag)
""")

# Isotropic material
E=1.0;
# nu=0.3;
nu=0.499999;
P=1.0;
L=16.0;
c=2.0;
h=1.0;
smult=0;
E1=E/(1-nu^2); nu1=nu/(1-nu); I=(2*c)^3/12;
uexx(x,y) = (P/(6*E1*I)/2*(-y).*(3*(L^2-(L-x).^2)+(2+nu1)*(y.^2-c^2)));
uexy(x,y) = (P/(6*E1*I)/2*(((L-x).^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x).*y.^2));

CTE = 0.0

tolerance = 0.00001*c

convergencestudy = FDataDict[]
for n = [1 2 4 8] #

    nL = 3*n # number of elements lengthwise
    nc = 2*n # number of elements through the wwith
    nh = n # number of elements through the thickness
    xs = collect(linspace(0.0, L, nL+1))
    ys = collect(linspace(0.0, h, nh+1))
    zs = collect(linspace(-c, +c, nc+1))
    fens,fes = T10blockx(xs, ys, zs)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
    # 0 cross-section surface  for the reactions
    sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR,  0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = eye(3, 3)

    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copy!(csmatout, csmat)
    end

    gr = SimplexRule(3, 4)

    region = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(fes, gr), material))

    lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    # println("lx0 = $(lx0)")
    ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
    ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
    lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
    lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
    # println("vcat(lx1, lx2) = $(vcat(lx1, lx2))")
    ex04 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>vcat(lx1, lx2) )
    ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
    # println("vcat(ly1, ly2) = $(vcat(ly1, ly2))")
    ey01 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly1, ly2) )

    function getfrcL!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        W=2/3*c^3;
        copy!(forceout, [0.0; 0.0; P/2/W*(c^2-XYZ[3]^2)])
    end

    function getfrc0!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        W=2/3*c^3;
        copy!(forceout, [P*L/W*XYZ[3]; 0.0; -P/2/W*(c^2-XYZ[3]^2)])
    end

    Trac0 = FDataDict("traction_vector"=>getfrc0!,
    "femm"=>FEMMBase(IntegData(subset(bfes, sshear0), SimplexRule(2, 3))))
    TracL = FDataDict("traction_vector"=>getfrcL!,
    "femm"=>FEMMBase(IntegData(subset(bfes, sshearL), SimplexRule(2, 3))))

    modeldata = FDataDict("fens"=>fens,
    "regions"=>[region],
    "essential_bcs"=>[ex01, ex02, ex03, ex04, ey01],
    "traction_bcs"=>[Trac0, TracL],
    "temperature_change"=>FDataDict("temperature"=>0.0)
    )
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

    u = modeldata["u"]
    geom = modeldata["geom"]

    Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
    utip = mean(u.values[Tipl, 3])
    println("Deflection: $(utip)")

    modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>[5])
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)

    modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>collect(1:6))
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    push!(convergencestudy, FDataDict(
        "elementsize"=> 1.0 / n,
        "fens"=>fens,
        "fes"=>fes,
        "geom"=>geom,
        "u"=>u,
        "femm"=>region["femm"],
        "integrationrule"=>region["femm"].integdata.integration_rule,
        "stressfields"=>stressfields,
        "tolerance"=>tolerance)
        )
end

File = "hughes_cantilever_stresses_incompressible_$(elementtag)"
open(File * ".jls", "w") do file
    serialize(file, convergencestudy)
end

ComputeErrorsModule.process(File)

true
