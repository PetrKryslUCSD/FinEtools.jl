module Planar_truss_examples
using FinEtools
using FinEtools.FENodeSetModule
using FinEtools.MeshExportModule
using LinearAlgebra

function Planar_truss()
    # Planar truss  structure loaded with concentrated forces  at some nodes.
    # The nodal displacements should be
    # 0         0
    # 0         0
    # 0.0213    0.0408
    # -0.0160    0.0462
    # 0.0427    0.1501
    # -0.0053    0.1661
    Area = 1.5
    E = 1.0e7 # Young's modulus
    nu = 0.0
    alpha = 0.0
    fens = FENodeSetModule.FENodeSet([0.0 0;
    0 40;
    40 0;
    40 40;
    80 0;
    80 40] )
    fes = FESetL2([1     3
    1     4
    2     4
    3     4
    3     5
    5     4
    6     4
    5     6])
    
    
    MR = DeforModelRed1D
    material = MatDeforElastIso(MR,  0.0, E, nu, alpha)
    # display(material )
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field
    setebc!(u, 1)
    setebc!(u, 2)
    applyebc!(u)
    numberdofs!(u)
    display(u)
    
    integdata = IntegData(fes, GaussRule(1, 1), (loc, conn, N) -> Area, false)
    femm = FEMMDeforLinear(MR, integdata, CSys(2, 1), material)
    display(femm.mcsys)
    K = stiffness(femm,  geom,  u)
    
    fi = ForceIntensity(vec([0 -2000.0]))
    lfemm = FEMMBase(IntegData(FESetP1(reshape([3], 1,1)), PointRule()))
    F = distribloads(lfemm,  geom,  u,  fi,  3);
    fi = ForceIntensity(vec([+2000.0 0]))
    lfemm = FEMMBase(IntegData(FESetP1(reshape([5], 1,1)), PointRule()))
    F = F + distribloads(lfemm,  geom,  u,  fi,  3);
    fi = ForceIntensity(vec([+4000.0 +6000.0]))
    lfemm = FEMMBase(IntegData(FESetP1(reshape([6], 1,1)), PointRule()))
    F = F + distribloads(lfemm,  geom,  u,  fi,  3);
    
    K = cholesky(K)
    U=  K\F
    scattersysvec!(u, U[:])
    
    sfld =  elemfieldfromintegpoints(femm, geom, u, :Cauchy, 1)
    display(sfld)
    println("Cauchy = $(sfld.values)")
    vfld =  elemfieldfromintegpoints(femm, geom, u, :vm, 1)
    display(vfld)
    
    File = "Planar_truss.vtk"
    MeshExportModule.vtkexportmesh(File, fens, fes;
    scalars=[("sx", sfld.values), ("vm", vfld.values)])
    @async run(`"paraview.exe" $File`)
    # try rm(File) catch end
    
end # Planar_truss

function allrun()
    println("#####################################################") 
    println("# Planar_truss ")
    Planar_truss()
    return true
end # function allrun

end # module Planar_truss_examples
