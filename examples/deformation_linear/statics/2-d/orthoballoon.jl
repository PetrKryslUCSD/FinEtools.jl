
module moorthobballoon

using FinEtools

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
# Note that the FinEtools objects needs to be created with the proper
# model-dimension reduction at hand.  In this case that is the axial symmetry
# assumption.
MR = DeforModelRed2DAxisymm


# Create the mesh and initialize the geometry.  First we are going
# to construct the block of elements with the first coordinate
# corresponding to the thickness in the radial direction, and the second
# coordinate is the thickness in the axial direction.

fens,fes = Q4block(rex-rin,pi/2,5,20);
bdryfes = meshboundary(fes);

icl = selectelem(fens, bdryfes, box=[0.,0.,0.,pi/2],inflate=tolerance);
for i=1:count(fens)
    r=rin+fens.xyz[i,1]; a=fens.xyz[i,2];
    fens.xyz[i,:]=[r*cos(a) r*sin(a)];
end


# now we create the geometry and displacement fields
geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

# the symmetry plane
l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
setebc!(u,l1,true, 2, 0.0)
# the axis of symmetry
l1 =selectnode(fens; box=[0 0 0 rex], inflate = tolerance)
setebc!(u,l1,true, 1, 0.0)

applyebc!(u)
numberdofs!(u)

# The traction boundary condition is applied in the radial
# direction.

el1femm =  FEMMBase(GeoD(subset(bdryfes,icl), GaussRule(1, 3), true))
function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copy!(forceout, XYZ/norm(XYZ)*p)
  return forceout
end
fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
F2= distribloads(el1femm, geom, u, fi, 2);

# Property and material
material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(2, 2), true), material)

K =stiffness(femm, geom, u)
U=  K\(F2)
scattersysvec!(u,U[:])

# Produce a plot of the radial stress component in the cylindrical
# coordinate system. Note that this is the usual representation of
# stress using nodal stress field.

fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 3)


File =  "orthoballoon_sigmaz.vtk"
vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
              vectors=[("u", u.values)])
@async run(`"paraview.exe" $File`)

show(minimum(fld.values))
show(maximum(fld.values))

end
#pub_thick_pipe_axi()
