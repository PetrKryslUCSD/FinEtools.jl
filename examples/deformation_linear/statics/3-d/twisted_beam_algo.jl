using FinEtools
using FinEtools.AlgoDeforLinearModule

println("""
The initially twisted cantilever beam is one of the standard test
problems for verifying the finite-element accuracy [1]. The beam is
clamped at one end and loaded either with unit in-plane or
unit out-of-plane force at the other. The centroidal axis of the beam is
straight at the undeformed  configuration, while its cross-sections are
twisted about the centroidal axis from 0 at the clamped end to pi/2 at
the free end.

Reference:
Zupan D, Saje M (2004) On "A proposed standard set of problems to test
finite element accuracy": the twisted beam. Finite Elements in Analysis
and Design 40: 1445-1451.
""")
function  Twisted_beam()
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 7;
  p =   1/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H20block(L,W,t, nl*ref,nw*ref,nt*ref)

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
  el1femm  = FEMMBase(GeoD(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR, GeoD(fes, GaussRule(3,2)),
            material))

  # Make model data
  modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1],
  "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
  geom = modeldata["geom"]
  u = modeldata["u"]

  # Extract the solution
  nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
  theutip = mean(u.values[nl,:],1)
  println("displacement  = $(theutip[dir]) as compared to converged $uex")
    println("normalized displacement  = $(theutip[dir]/uex*100) %")

  # Write out mesh with displacements
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam")
  modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :Cauchy, "component"=> :xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)

  # Write out mesh with von Mises stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  vm  = modeldata["postprocessing"]["exported_fields"][1]
  println("extremes of von Mises: $([minimum(vm.values),   maximum(vm.values)])")

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)

  println("Done")
  true

end
Twisted_beam()
