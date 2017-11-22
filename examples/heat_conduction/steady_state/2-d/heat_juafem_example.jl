using FinEtools
using FinEtools.MeshExportModule
using FinEtools.IntegDataModule: integrationdata
using FinEtools.FESetModule.gradN!

println("""Heat conduction example from JuAFEM.""")
t0 = time()

A = 2.0
thermal_conductivity = eye(2,2); # conductivity matrix
function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  forceout[1] = 1.0; #heat source
end
N = 1000;

println("Mesh generation")
@time fens,fes = Q4block(A, A, N, N)
fens.xyz[:,1] .-= A/2
fens.xyz[:,2] .-= A/2

geom = NodalField(fens.xyz)
Temp = NodalField(zeros(size(fens.xyz,1),1))

println("Searching nodes  for BC")
@time l1 = selectnode(fens; box=[-A/2 -A/2 -A/2 A/2], inflate = 1.0/N/100.0)
@time l2 = selectnode(fens; box=[A/2 A/2 -A/2 A/2], inflate = 1.0/N/100.0)
@time l3 = selectnode(fens; box=[-A/2 A/2 -A/2 -A/2], inflate = 1.0/N/100.0)
@time l4 = selectnode(fens; box=[-A/2 A/2 A/2 A/2], inflate = 1.0/N/100.0)
List = vcat(l1, l2, l3, l4);
setebc!(Temp, List, true, 1, 0.0)
applyebc!(Temp)

numberdofs!(Temp)

t1 = time()

m = MatHeatDiff(thermal_conductivity)
femm = FEMMHeatDiff(IntegData(fes, GaussRule(2, 2)), m)

println("Conductivity")
@time K=conductivity(femm, geom, Temp)

println("Internal heat generation")
fi = ForceIntensity(FFlt, 1, getsource!);
@time F1 = distribloads(femm, geom, Temp, fi, 3);

println("Factorization")
@time K = cholfact(K)
println("Solution of the factorized system")
@time U =  K\(F1)
scattersysvec!(Temp, U[:])


println("Total time elapsed = $(time() - t0) [s]")
println("Solution time elapsed = $(time() - t1) [s]")

println("Maximum temperature = $(maximum(Temp.values)) ")
# using MeshExportModule

# File =  "a.vtk"
# MeshExportModule.vtkexportmesh(File, fes.conn, hcat(geom.values,  Temp.values), MeshExportModule.Q4; scalars=[("Temperature", Temp.values)])
# @async run(`"paraview.exe" $File`)

function  buffers(self::FEMMHeatDiff, geom::NodalField{FFlt}, temp::NodalField{FFlt})
    # Constants
    nfes = count(self.integdata.fes); # number of finite elements in the set
    ndn = ndofs(temp); # number of degrees of freedom per node
    nne = nodesperelem(self.integdata.fes); # number of nodes for element
    sdim = ndofs(geom);   # number of space dimensions
    mdim = manifdim(self.integdata.fes); # manifold dimension of the element
    Kedim = ndn*nne;      # dimension of the element matrix
    elmat = zeros(FFlt, Kedim, Kedim); # buffer
    elvec = zeros(FFlt, Kedim); # buffer
    elvecfix = zeros(FFlt, Kedim); # buffer
    conn = zeros(FInt, nne); # buffer
    x = zeros(FFlt, nne, sdim); # buffer
    dofnums = zeros(FInt, 1, Kedim); # buffer
    loc = zeros(FFlt, 1, sdim); # buffer
    J = eye(FFlt, sdim, mdim); # buffer
    RmTJ = zeros(FFlt, mdim, mdim); # buffer
    gradN = zeros(FFlt, nne, mdim); # buffer
    kappa_bargradNT = zeros(FFlt, mdim, nne); # buffer
    return conn, x, dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat, elvec, elvecfix
end
function justfun(self::FEMMHeatDiff, geom::NodalField{FFlt},  temp::NodalField{FFlt}) 
    integdata = self.integdata
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(integdata);
    # Prepare assembler and buffers
    conn, x, dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat = buffers(self, geom, temp)
    for i = 1:count(integdata.fes) # Loop over elements
        getconn!(integdata.fes, conn, i);
        gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
        for j=1:npts # Loop over quadrature points
            At_mul_B!(loc, Ns[j], x);# Quadrature point location
            At_mul_B!(J, x, gradNparams[j]); # Jacobian matrix
            Jac = Jacobianvolume(integdata, J, loc, conn, Ns[j]);
            # updatecsmat!(self.mcsys, loc, J, integdata.fes.label[i]);
            # At_mul_B!(RmTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
            # gradN!(integdata.fes, gradN, gradNparams[j], RmTJ);
            gradN!(integdata.fes, gradN, gradNparams[j], J);
        end # Loop over quadrature points
    end # Loop over elements
    return nothing;
end

println("Just fun")
@time justfun(femm, geom, Temp)

true
