
module mmLE1NAFEMSstressESNICET4
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Radius = 1.0*phun("m")
    Thickness = 0.1*phun("m")
    n = 8; # number of elements per side
    tolerance = 1.0/n/1000.; # Geometrical tolerance

    fens,fes = T4block(1.0, pi/2, Thickness, n, n*2, 3)

    bdryfes = meshboundary(fes);
    icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
    for i=1:count(fens)
        t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
        fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
    end

    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 1, 0.0)
    l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
    setebc!(u,l1,true, 3, 0.0)

    applyebc!(u)
    numberdofs!(u)

    el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), SimplexRule(2, 1)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
        forceout .=    vec(p*pt/norm(pt));
        return forceout
    end
    fi = ForceIntensity(FFlt, 3, pfun);
    F2 = distribloads(el1femm, geom, u, fi, 2);

    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D

    material = MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)

    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)

    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K\(F2)
    scattersysvec!(u, U[:])

    # File =  "a.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.T4; vectors = [("u", u.values)])
    # @async run(`"paraview.exe" $File`)

    nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, 0.0], inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.10215,0] [MM]")
    @test norm(thecorneru - [-0.115443 0.0 0.0]) < 1.0e-4


    fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging)
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    # File =  "a.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.T4; vectors = [("u", u.values)], scalars = [("sig", fld.values)])
    # @async run(`"paraview.exe" $File`)
    @test abs(fld.values[nl,1][1]/phun("MPa") - 91.07617647479451) < 1.0e-3

end
end
using .mmLE1NAFEMSstressESNICET4
mmLE1NAFEMSstressESNICET4.test()


