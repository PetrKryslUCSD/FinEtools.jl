module m # pub_thick_pipe_axi()

using JFinEALE

# NAFEMS LE11 benchmark with Q4 elements.
# # This is a test recommended by the National Agency for Finite Element
# # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
# # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.  
# #
# # Target solution: Direct stress,  = –105 MPa at point A.
#function  LE11NAFEMS()
# Parameters:
Ea= 210000*phun("MEGA*Pa")
nua= 0.3;
alphaa=2.3e-4;              # thermal expansion coefficient
sigmaA=-105*phun("MEGA*Pa")
nref= 1;                        # how many times should we refine the mesh?
X=[1.     0.;#A
   1.4    0.;#B
   0.995184726672197   0.098017140329561;
   1.393258617341076 0.137223996461385;
   0.980785  0.195090;#
   1.37309939 0.27312645;
   0.956940335732209   0.290284677254462
   1.339716470025092 0.406398548156247
   0.9238795  0.38268;#C
   1.2124  0.7;#D
   0.7071  0.7071;#E
   1.1062  1.045;#F
   0.7071  (0.7071+1.79)/2;#(E+H)/2
   1.      1.39;#G
   0.7071  1.79;#H
   1.      1.79;#I
   ]*phun("M")
   tolerance =1.e-6*phun("M")
   ##
   # Note that the material object needs to be created with the proper
   # model-dimension reduction in mind.  In this case that is the axial symmetry
   # assumption.
   mr=DeformationModelReduction2DAxisymm



fens=FENodeSet(X);
fes=FESetQ4(conn=[1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
for ref=1:nref
    fens,fes=Q4refine(fens,fes);
    list=selectnode(fens,distance=1.0+0.1/2^nref, from=[0. 0.], inflate=tolerance);
    fens.xyz[list,:]= JFinEALE.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
end 
fes.axisymm=(mr== DeformationModelReduction2DAxisymm)  # note that this reflects the chosen model reduction

#     File =  "mesh.vtk"
# vtkexportmesh(File, fens, fes)

# now we create the geometry and displacement fields
geom = NodalField(name ="geom",data =fens.xyz)
u = NodalField(name ="u",data =zeros(size(fens.xyz,1),2)) # displacement field

# Apply EBC's
l1=selectnode(fens,box=[-Inf Inf 0 0],inflate=tolerance)
setebc!(u,l1,trues(length(l1)),l1*0+2,[0.0])    
l1=selectnode(fens,box=[-Inf Inf 1.79  1.79],inflate=tolerance)
setebc!(u,l1,trues(length(l1)),l1*0+2,[0.0])    
applyebc!(u)
numberdofs!(u)

# Temperature field
dT =NodalField(name="dT",data=reshape(fens.xyz[:,1]+fens.xyz[:,2],size(fens.xyz,1),1));


# Property and material
material=MaterialDeformationLinear (PropertyDeformationLinearIso(E=Ea,nu=nua,CTE=alphaa))

femm = FEMMDeformationLinear(FEMMBase(fes, GaussRule(order=2,dim=2)), material)

K =stiffness(mr, femm, geom, u)
F = thermalstrainloads(mr, femm, geom, u, dT)
#K=cholfact(K)
U=  K\F
scattersysvec!(u,U[:])

nA =selectnode(fens,box=JFFlt[1.0  1.0 0.0 0.0],inflate=tolerance);

fld= fieldfromintegpoints(mr, femm, geom, u, dT, :Cauchy, 2)


File =  "LE11NAFEMS_Q4_sigmay.vtk"
vtkexportmesh(File, fens, fes; scalars=fld.values,scalars_name ="sigmay", vectors=u.values,vectors_name="u")

sA = fld.values[nA]/phun("MEGA*Pa")
sAn = fld.values[nA]/sigmaA
println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")

fen2fe =FENodeToFEMap(fes.conn,nnodes(geom))
function inspector(idat,out,loc,pc)
    println("loc=$(  loc  )")
    println(" : $(  out'  )")
    return idat
end

inspectintegpoints(mr,femm, geom, u, dT,  fen2fe.map[nA[1]],
                   inspector, []; output=:Cauchy)

#finealemesh(fens,fes,"meshmfile")

# end
# LE11NAFEMS()

end

