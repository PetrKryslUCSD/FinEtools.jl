module scratch

using JFinEALE


## Harmonic forced vibration analysis of simply-supported thin  (solid) plate
#

##
# Link to the  <matlab:edit('pub_TEST13H_vibration') m-file>.
#

## Description
#
# Harmonic forced vibration problem is solved for a homogeneous square plate,
# simply-supported on the circumference.
# This is the TEST 13H from the Abaqus v 6.12 Benchmarks manual.
# The test is recommended by the National Agency for Finite Element Methods and Standards (U.K.):
# Test 13 from NAFEMS “Selected Benchmarks for Forced Vibration,” R0016, March 1993.

##
# The plate is discretized with hexahedral solid elements. The simple support
# condition is approximated by distributed rollers on the boundary.
# Because only the out of plane displacements are prevented, the structure
# has three rigid body modes in the plane of the plate.

##
# The nonzero benchmark frequencies are (in hertz): 2.377, 5.961, 5.961,
# 9.483, 12.133, 12.133, 15.468, 15.468 [Hz]. 

println("""
        Harmonic forced vibration analysis of simply-supported thin  (solid) plate.
        The nonzero benchmark frequencies are (in hertz):
        2.377, 5.961, 5.961, 9.483, 12.133, 12.133, 15.468, 15.468 [Hz]. 
        """) 
        t0 = time()


##
# Define the material properties.

# Parameters:
E = 200*phun("GIGA*PA");
nu = 0.3;
rho= 8000*phun("KG/M^3");

##
# Geometrical dimensions of the plate (the full structure).
L =10*phun("M");# span of the plate
t =0.05*phun("M");# thickness of the plate

##
# Applied pressure on the top surface
Q=100.*phun("N/M^2")

##
# The chosen mesh parameters. This is the  coarse mesh as specified in the benchmark.
nL= 8;# number of elements span wise
nt = 1;# number of elements through the thickness

neigvs=20                   # how many eigenvalues
OmegaShift=(0.01*2*pi)^2;

##
# The mesh is generated. The chosen elements are the serendipity hexahedra.
# The model is 3-D.
Modelreduction=DeformationModelReduction3D
fens,fes =H20block(L,L,t,nL,nL,nt);;

geom = NodalField(name ="geom",data =fens.xyz)
u = NodalField(name ="u",data =zeros(JFCplxFlt,size(fens.xyz,1),3)) # displacement field


##
#
# The support conditions approximate simply-supported edges.  All the
# sides of the plate are fixed in the transverse direction  (Z
# displacement).
nodelist = selectnode(fens, box= [0.0,0.0,-Inf,Inf,-Inf,Inf], inflate=0.001*t);
for i in nodelist
    setebc!(u,i, true, 3, 0.0)
end
nodelist = selectnode(fens, box= [L,L,-Inf,Inf,-Inf,Inf], inflate=0.001*t);
for i in nodelist
    setebc!(u,i, true, 3, 0.0)
end                                                                                
nodelist = selectnode(fens, box= [0.0,L,0.0,0.0,-Inf,Inf], inflate=0.001*t)
for i in nodelist
    setebc!(u,i, true, 3, 0.0)
end
nodelist = selectnode(fens, box= [0.0,L,L,L,-Inf,Inf], inflate=0.001*t)
for i in nodelist
    setebc!(u,i, true, 3, 0.0)
end                                                                                                                                                                                                                                                                                                        

numberdofs!(u)

p=PropertyDeformationLinearIso(rho,E,nu)
material=MaterialDeformationLinear (p)

# Stiffness matrix uses uniformly reduced integration
femm = FEMMDeformationLinear(FEMMBase(fes, GaussRule(order=2,dim=3)), material)
K =stiffness(Modelreduction, femm, geom, u)
# Mass matrix  is based on full integration
femm = FEMMDeformationLinear(FEMMBase(fes, GaussRule(order=3,dim=3)), material)
M =mass(Modelreduction, femm, geom, u)

# Solve the eigenvalue problem
d,v,nev,nconv =eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
d = d - OmegaShift;
fs=real(sqrt(complex(d)))/(2*pi)
println("Eigenvalues: $fs [Hz]")

mode=7
scattersysvec!(u,v[:,mode])
File =  "TEST13H_vibration.vtk"
vtkexportmesh (File, fens, fes; vectors=u.values,vectors_name ="mode$mode")


## 
# Now we are going to modify the above three-vibration model to incorporate
# harmonic forcing.  The entire  top surface of the plate  is loaded with
# uniform pressure.
##
# Define the traction load of  100 N/m^2  over the entire surface.  This traction load has harmonic
# dependence but its distribution does not change as a function of
# frequency. 
bfes= meshboundary(fes);
topl =selectelem (fens,bfes,box= [-Inf,Inf,-Inf,Inf,t,t], inflate=t/100);
el1femm =  FEMMBase(subset(bfes,topl), GaussRule(order=2,dim=2))
fi = ForceIntensity(complex([0.0,0.0,Q]));
F= distribloads(el1femm, geom, u, fi, 2);

##
# Compute the parameters of Rayleigh damping. For the two selected
# frequencies we have the relationship between the damping ratio and
# the Rayleigh parameters
##
# $\xi_m=a_0/\omega_m+a_1\omega_m$

##
# where $m=1,2$.  Solving for the Rayleigh parameters $a_0,a_1$ yields:
zeta1= 0.02; zeta2  =0.02;
o1 =2*pi*2.377;  o2 =2*pi*15;
Rayleigh_mass = 2*(o1*o2)/(o2^2-o1^2)*(o2*zeta1-o1*zeta2);# a0
Rayleigh_stiffness = 2*(o1*o2)/(o2^2-o1^2)*(-1/o2*zeta1+1/o1*zeta2);# a1
# Damping matrix
C=Rayleigh_mass*M+Rayleigh_stiffness*K

##
# These are the frequencies at which to evaluate the frequency response
# function. Note that we are taking one of the points as the predicted fundamental frequency.
frequencies = [linspace(0,2.377,15),linspace(2.377,15,30)];

##
# The function below  will be called with each computed displacement
# from within the solver.  The amplitude of the deflection at the
# midpoint in the direction of the load will be saved  for each
# frequency.
midpoint=selectnode (fens, box=[L/2 L/2 L/2 L/2 0 0],inflate=t/100);
midpointuz=zeros(JFCplxFlt,length(frequencies))


for i=1:length(frequencies)
    f=frequencies[i]
    omega=2*pi*f
    A = lufact(K+1.0im*omega*C-omega^2*M)
    U = A\F
    scattersysvec!(u,U)
    midpointuz[i]=u.values[midpoint,3][1]
    #println("$( u.values[midpoint,:]  )")
end

using Winston

using Winston

pl = FramedPlot(title="Thin plate midpoint FRF",xlabel="Frequency [Hz]",ylabel="Midpoint  displacement amplitude [mm]")
setattr(pl.frame, draw_grid=true)
add(pl, Curve(frequencies,real(midpointuz)/phun("mm"), color="blue"))
add(pl, Curve(frequencies,imag(midpointuz)/phun("mm"), color="red"))
display(pl)

pl = FramedPlot(title="Thin plate midpoint FRF",xlabel="Frequency [Hz]",ylabel="Phase shift of the midpoint FRF")
setattr(pl.frame, draw_grid=true)
add(pl, Curve(frequencies,atan2(imag(midpointuz),real(midpointuz) )/pi*180, color="green"))
display(pl)



end
