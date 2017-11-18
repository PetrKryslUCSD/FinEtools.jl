
module mocylpull

using FinEtools

# Cylinder  pulled by enforced displacement, axially symmetric model


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
Length = 1*rex
ua = -0.05*Length
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

fens,fes = Q4block(rex-rin,Length,5,20);
fens.xyz[:, 1] += rin
bdryfes = meshboundary(fes);

# now we create the geometry and displacement fields
geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

# the symmetry plane
l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
setebc!(u,l1,true, 2, 0.0)
# The other end
l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
setebc!(u,l1,true, 2, ua)

applyebc!(u)
numberdofs!(u)
println("Number of degrees of freedom = $(u.nfreedofs)")

# Property and material
material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)

K =stiffness(femm, geom, u)
F = nzebcloadsstiffness(femm, geom, u)
U=  K\(F)
scattersysvec!(u,U[:])

# Produce a plot of the radial stress component in the cylindrical
# coordinate system. Note that this is the usual representation of
# stress using nodal stress field.

fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2)
println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")

File =  "orthoballoon_sigmaz.vtk"
vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
              vectors=[("u", u.values)])
@async run(`"paraview.exe" $File`)


end
#pub_thick_pipe_axi()
