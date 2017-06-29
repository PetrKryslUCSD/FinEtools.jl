using JFinEALE

println("""
% Vibration modes of unit cube  of almost incompressible material.
%
% Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
% tetrahedral. International Journal for Numerical Methods in
% Engineering 67: 841-867.""") 
t0 = time()


E = 1*phun("PA");
nu = 0.499;
rho= 1*phun("KG/M^3");
a=1*phun("M"); b=a; h= a;
n1=10;# How many element edges per side?
na= n1; nb= n1; nh =n1;
neigvs=20                   # how many eigenvalues
omega_shift=(0.01*2*pi)^2;

Modelreduction=DeformationModelReduction3D
fens,fes =H20block(a,b,h, na,nb,nh)

# Make the region
p=PropertyDeformationLinearIso(rho,E,nu)
region1=dmake(fes=fes, property=p,
              integration_rule_stiffness=GaussRule(order=2, dim=3),
              integration_rule_mass=GaussRule(order=3, dim=3))

# Make model data
modeldata= dmake(modelreduction=DeformationModelReduction3D,
                 omega_shift=omega_shift,
                 neigvs=neigvs,
                 fens= fens,region=[region1],
                 boundary_conditions=dmake() # the cube is free-floating
                 );

# Solve
modeldata=JFinEALE.DeformationLinearAlgorithmModule.modal(modeldata)

fs=modeldata["omega"]/(2*pi)
println("Eigenvalues: $fs [Hz]")

modeldata["postprocessing"] = dmake(file="unit_cube_mode1",mode=10)
modeldata=JFinEALE.DeformationLinearAlgorithmModule.exportmode(modeldata)

true

    
