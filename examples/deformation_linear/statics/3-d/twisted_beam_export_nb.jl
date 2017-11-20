using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule

E = 0.29e8;
nu = 0.22;
W = 1.1;
L = 12.;
t =  0.32;
nl = 2; nt = 1; nw = 1; ref = 4;
p =   1/W/t;
#  Loading in the Z direction
loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
#   Loading in the Y direction
#loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
tolerance  = t/1000;

fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)

# Reshape into a twisted beam shape
for i = 1:count(fens)
  a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
  fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
end

# Clamped end of the beam
l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

# Traction on the opposite edge
boundaryfes  =   meshboundary(fes);
Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
el1femm  = FEMMBase(IntegData(subset(boundaryfes,Toplist), GaussRule(2, 2)))
flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


# Make the region
MR = DeforModelRed3D
material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3,2)),
material))

# Make model data
modeldata =  FDataDict(
"fens"=> fens, "regions"=>  [region1],
"essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])


AE = AbaqusExporter("twisted_beam");
# AE.ios = STDOUT;
HEADING(AE, "Twisted beam example");
PART(AE, "part1");
END_PART(AE);
ASSEMBLY(AE, "ASSEM1");
INSTANCE(AE, "INSTNC1", "PART1");
NODE(AE, fens.xyz);
ELEMENT(AE, "c3d8rh", "AllElements", 1, region1["femm"].integdata.fes.conn)
ELEMENT(AE, "SFM3D4", "TractionElements",
1+count(region1["femm"].integdata.fes), flux1["femm"].integdata.fes.conn)
NSET_NSET(AE, "l1", l1)
ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl");
SURFACE_SECTION(AE, "TractionElements")
END_INSTANCE(AE);
END_ASSEMBLY(AE);
MATERIAL(AE, "elasticity")
ELASTIC(AE, E, nu)
SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
STEP_PERTURBATION_STATIC(AE)
BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 1)
BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 2)
BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 3)
DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec(flux1["traction_vector"]))
END_STEP(AE)
close(AE)
