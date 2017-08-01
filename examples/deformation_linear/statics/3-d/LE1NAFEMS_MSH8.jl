using FinEtools
using FinEtools.MeshExportModule

println("LE1NAFEMS, 3D version."        )
t0 = time()

E = 210e3*phun("MEGA*PA");# 210e3 MPa
nu = 0.3;
p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
sigma_yD= 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
Radius= 1.0*phun("m")
Thickness = 0.1*phun("m")
n=32; # number of elements per side
tolerance=1.0/n/1000.;#Geometrical tolerance

fens,fes =Q4block(1.0,pi/2, n, n*2)

bdryfes = meshboundary(fes);
icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness],inflate=tolerance);
for i=1:count(fens)
    t=fens.xyz[i,1]; a=fens.xyz[i,2];
    fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a) (t*2.75+(1-t)*1)*sin(a)];
end


geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

l1 =selectnode(fens; box=[0.0 Inf 0.0 0.0], inflate = tolerance)
setebc!(u,l1,true, 2, 0.0)
l1 =selectnode(fens; box=[0.0 0.0 0.0 Inf], inflate = tolerance)
setebc!(u,l1,true, 1, 0.0)

applyebc!(u)
numberdofs!(u)


el1femm =  FEMMBase(GeoD(subset(bdryfes,icl), GaussRule(2, 2)))
function pfun(x::JFFltMat,J::JFFltMat,l::JFInt)
    pt= [2.75/3.25*x[1], 3.25/2.75*x[2], 0.0]
    return   (p*pt'/norm(pt));
end
fi = ForceIntensity(JFFlt,pfun);
F2= distribloads(el1femm, geom, u, fi, 2);


p=PropertyDeformationLinearIso(E,nu)
material=MaterialDeformationLinear (p)

femm = FEMMDeformationLinear(FEMMBase(fes, GaussRule(order=2,dim=2)), material)

# The geometry field now needs to be associated with the FEMM
femm = associategeometry!(femm, geom)

K =stiffness(femm, geom, u)
K=cholfact(K)
U=  K\(F2)
scattersysvec!(u,U[:])

nl=selectnode (fens, box=[2.0,2.0,0.0,0.0],inflate=tolerance);
thecorneru=zeros(JFFlt,1,2)
gathervaluesasmat!(u,thecorneru,nl);
thecorneru=thecorneru/phun("mm")
println("$(time()-t0) [s];  displacement =$(thecorneru) [MM] as compared to reference [-0.10215,0] [MM]")



File =  "a.vtk"
vtkexportmesh (File, fes.conn, geom.values,
               FinEtools.MeshExportModule.Q4; vectors=u.values,vectors_name ="u")

true
