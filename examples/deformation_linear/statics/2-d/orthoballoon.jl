
module m # pub_thick_pipe_axi()

using JFinEALE

# Orthotropic balloon inflation, axially symmetric model


# Parameters:
E1=1.0;
E2=1.0;
E3=3.0;
nu12=0.29;
nu13=0.29;
nu23=0.19;
G12=0.3;
G13=0.3;
G23=0.3;
p= 0.15;
rin=1.;
rex =1.2;
tolerance=rin/1000.

##
# Note that the JFinEALE objects needs to be created with the proper
# model-dimension reduction at hand.  In this case that is the axial symmetry
# assumption.
mr=DeformationModelReduction2DAxisymm


# Create the mesh and initialize the geometry.  First we are going
# to construct the block of elements with the first coordinate
# corresponding to the thickness in the radial direction, and the second
# coordinate is the thickness in the axial direction.

fens,fes=Q4block(rex-rin,pi/2,5,20);
fes.axisymm=(mr== DeformationModelReduction2DAxisymm)           # note that this reflects the chosen model reduction
bdryfes = meshboundary(fes);
bdryfes.axisymm=(mr== DeformationModelReduction2DAxisymm)

icl = selectelem(fens, bdryfes, box=[0.,0.,0.,pi/2],inflate=tolerance);
for i=1:count(fens)
    r=rin+fens.xyz[i,1]; a=fens.xyz[i,2];
    fens.xyz[i,:]=[r*cos(a) r*sin(a)];
end


# now we create the geometry and displacement fields
geom = NodalField(name ="geom",data =fens.xyz)
u = NodalField(name ="u",data =zeros(size(fens.xyz,1),2)) # displacement field

# the symmetry plane
l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
setebc!(u,l1,trues(length(l1)),l1*0+2,[0.0])
# the axis of symmetry
l1 =selectnode(fens; box=[0 0 0 rex], inflate = tolerance)
setebc!(u,l1,trues(length(l1)),l1*0+1,[0.0])

applyebc!(u)
numberdofs!(u)

# The traction boundary condition is applied in the radial
# direction.
    
el1femm =  FEMMBase(subset(bdryfes,icl), GaussRule(order=3,dim=1))
fi = ForceIntensity(JFFlt,(x,J,l)->p*x/norm(x));
F2= distribloads(el1femm, geom, u, fi, 2);

# Property and material
material=MaterialDeformationLinear (PropertyDeformationLinearOrtho(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23))

femm = FEMMDeformationLinear(FEMMBase(fes, GaussRule(order=2,dim=2)), material)

K =stiffness(mr, femm, geom, u)
#K=cholfact(K)
U=  K\(F2)
scattersysvec!(u,U[:])

# Produce a plot of the radial stress component in the cylindrical
# coordinate system. Note that this is the usual representation of
# stress using nodal stress field.

fld= fieldfromintegpoints(mr, femm, geom, u, :Cauchy, 3)


File =  "orthoballoon_sigmaz.vtk"
vtkexportmesh(File, fens, fes; scalars=fld.values,scalars_name ="sigmaz",
              vectors=[u.values zeros(JFFlt,size(u.values,1))])

                              
            

end
#pub_thick_pipe_axi()
