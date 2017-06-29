using JFinEALE

println("LE1NAFEMS, plane stress, Q8 elements."        )        
t0 = time()

E = 210e3*phun("MEGA*PA");# 210e3 MPa
nu = 0.3;
p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
sigma_yD= 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
Radius= 1.0*phun("m")
n=3; # number of elements per side
tolerance=1.0/n/1000.;#Geometrical tolerance

fens,fes =Q8block(1.0,pi/2, n, n*2)
setotherdimension!(fes,1.0)

bdryfes = meshboundary(fes);
icl = selectelem(fens, bdryfes, box=[1.0,1.0,0.0,pi/2],inflate=tolerance);
for i=1:count(fens)
    t=fens.xyz[i,1]; a=fens.xyz[i,2];
    fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a) (t*2.75+(1-t)*1)*sin(a)];
end


geom = NodalField(name ="geom",data =fens.xyz)
u = NodalField(name ="u",data =zeros(size(fens.xyz,1),2)) # displacement field

l1 =selectnode(fens; box=[0.0 Inf 0.0 0.0], inflate = tolerance)
setebc!(u,l1,trues(length(l1)),l1*0+2,[0.0])
l1 =selectnode(fens; box=[0.0 0.0 0.0 Inf], inflate = tolerance)
setebc!(u,l1,trues(length(l1)),l1*0+1,[0.0])    

applyebc!(u)
numberdofs!(u)

el1femm =  FEMMBase(subset(bdryfes,icl), GaussRule(order=2,dim=1))
function pfun(x::JFFltMat,J::JFFltMat,l::JFInt)
    pt= [2.75/3.25*x[1],3.25/2.75*x[2]]
    return   (p*pt'/norm(pt));
end
fi = ForceIntensity(JFFlt,pfun);
F2= distribloads(el1femm, geom, u, fi, 2);


p=PropertyDeformationLinearIso(E,nu)
material=MaterialDeformationLinear (p)

femm = FEMMDeformationLinear(FEMMBase(fes, GaussRule(order=2,dim=2)), material)

K =stiffness(DeformationModelReduction2DStress, femm, geom, u)
K=cholfact(K)
U=  K\(F2)
scattersysvec!(u,U[:])

nl=selectnode (fens, box=[2.0,2.0,0.0,0.0],inflate=tolerance);
thecorneru=zeros(JFFlt,1,2)
gathervaluesasmat!(u,thecorneru,nl);
thecorneru=thecorneru/phun("mm")
println("$(time()-t0) [s];  displacement =$(thecorneru) [MM] as compared to reference [-0.10215,0] [MM]")

function inspector(idat,out,loc,pc)
    println("$(  loc  ): $(  out  )")
    return idat
end
felist=[1]
inspectintegpoints(DeformationModelReduction2DStress,femm, 
                   geom, u,
                   felist,
                   inspector, [])

fld= fieldfromintegpoints(DeformationModelReduction2DStress,femm, 
                          geom, u, :Cauchy, 2)
                          
                          using JFinEALE.MeshExportModule

File =  "a.vtk"
vtkexportmesh(File, fens, fes; scalars=fld.values,scalars_name ="sigmay")

true
