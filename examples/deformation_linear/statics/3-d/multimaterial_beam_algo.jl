using JFinEALE
using JFinEALE.AlgoDeformationLinearModule
using JFinEALE.MeshExportModule

println("""
Multi-material beam. Rubber-like and metal-like halfs, 
clamped, with shear traction at free end.
""")
function  multimaterial_beam()
    E1=0.29e3;
    nu1=0.49;
    E2=0.4e4;
    nu2=0.3;
    W=4.1;
    L=12.;
    t=6.5;
    nl=2; nt=1; nw=1; ref=9;
    p=  200./W/t;
    #  Loading in the Z direction
    loadv=[0;0;p];dir=3; 
    tolerance =t/1000;
    
    fens,fes =H20block(L,W,t, nl*ref,nw*ref,nt*ref)

    # Clamped end of the beam
    l1 =selectnode(fens; box=[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    essential1=dmake(node_list=l1,component=1:3,displacement=0.0)

    # Traction on the opposite edge
    boundaryfes =  meshboundary(fes);
    Toplist  =selectelem(fens,boundaryfes, box= [L L -Inf Inf -Inf Inf], inflate=  tolerance);
    el1femm =  FEMMBase(subset(boundaryfes,Toplist), GaussRule(order=2,dim=2))
    flux1=dmake(femm=el1femm,traction_vector=loadv)

    r1list  =selectelem(fens,fes, box= [0 L/2. -Inf Inf -Inf Inf], inflate=  tolerance);
    r2list  =selectelem(fens,fes, box= [L/2. L -Inf Inf -Inf Inf], inflate=  tolerance);
    
    # Make region 1
    p1=PropertyDeformationLinearIso(E1,nu1)
    region1=dmake(fes=subset(fes,r1list), property=p1, integration_rule=GaussRule(order=2,dim=3))

    # Make region 2
    p2=PropertyDeformationLinearIso(E2,nu2)
    region2=dmake(fes=subset(fes,r2list), property=p2, integration_rule=GaussRule(order=2,dim=3))

    # Model reduction type
    mr=DeformationModelReduction3D

    # Make model data
    modeldata= dmake(modelreduction=mr,
                     fens= fens,region=[region1 region2],
                     boundary_conditions=dmake(traction=[flux1],essential=[essential1]));

    # Call the solver
    modeldata=JFinEALE.AlgoDeformationLinearModule.linearstatics(modeldata)
    geom=modeldata["geom"]
    u=modeldata["u"]

    # Write out mesh with displacements
    dadd!(modeldata,postprocessing=dmake(file=  "multimaterial_beam"))
    modeldata=JFinEALE.AlgoDeformationLinearModule.exportdeformation(modeldata)

    # Write out mesh with stresses
    dadd!(modeldata,postprocessing=dmake(file=  "multimaterial_beam", output=:Cauchy, component=:xy))
    modeldata=JFinEALE.AlgoDeformationLinearModule.exportstress(modeldata)

    # Write out mesh with stresses
    dadd!(modeldata,postprocessing=dmake(file=  "multimaterial_beam", output=:Cauchy, component=:xz))
    modeldata=JFinEALE.AlgoDeformationLinearModule.exportstress(modeldata)

    # Write out mesh with von Mises stresses
    dadd!(modeldata,postprocessing=dmake(file=  "multimaterial_beam", output=:vm))
    modeldata=JFinEALE.AlgoDeformationLinearModule.exportstress(modeldata)

    # Write out mesh with von Mises stresses, elementwise
    dadd!(modeldata,postprocessing=dmake(file=  "multimaterial_beam-ew", output=:vm))
    modeldata=JFinEALE.AlgoDeformationLinearModule.exportstresselementwise(modeldata)

    # Write out mesh with von Mises stresses, elementwise
    dadd!(modeldata,postprocessing=dmake(file=  "multimaterial_beam-ew", output=:Cauchy, component=:xz))
    modeldata=JFinEALE.AlgoDeformationLinearModule.exportstresselementwise(modeldata)

    true

end
multimaterial_beam()
    
