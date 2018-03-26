module mmmmPoiss_06122017
using FinEtools
using Compat.Test
import LinearAlgebra: cholfact
function test()

  # println("""
  #
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit square, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
  # in a grid of 1000 x 1000 edges (2M triangles, 1M degrees of freedom).
  # Version: 05/29/2017
  # """
  # )
  t0 = time()

  A = 1.0 # dimension of the domain (length of the side of the square)
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
  Q = -6.0; # internal heat generation rate
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = Q; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 100;# number of subdivisions along the sides of the square domain


  # println("Mesh generation")
  fens,fes =T3block(A, A, N, N)

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  # println("Searching nodes  for BC")
  l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
  l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
  l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
  l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
  List = vcat(l1, l2, l3, l4)
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)
  numberdofs!(Temp)

  t1 = time()

  material = MatHeatDiff(thermal_conductivity)

  femm = FEMMHeatDiff(IntegData(fes, TriRule(1), 100.), material)


  # println("Conductivity")
  K = conductivity(femm, geom, Temp)
  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
  fi = ForceIntensity(FFlt[Q]);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")
  K = cholfact(K)
  # println("Solution of the factorized system")
  U = K\(F1+F2)
  scattersysvec!(Temp,U[:])

  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  Error= 0.0
  for k=1:size(fens.xyz,1)
    Error = Error.+abs.(Temp.values[k,1].-tempf(reshape(fens.xyz[k,:], (1,2))))
  end
  # println("Error =$Error")


  # File =  "a.vtk"
  # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

  @test Error[1]<1.e-5

  true
end
end
using .mmmmPoiss_06122017
mmmmPoiss_06122017.test()

module mmmmannulus_Q4_example_algo
using FinEtools
using Compat.Test
function test()
  # println("""
  # Annular region, ingoing and outgoing flux. Temperature at one node prescribed.
  # Minimum/maximum temperature ~(+/-)0.500.
  # Mesh of serendipity quadrilaterals.
  # This version uses the FinEtools algorithm module.
  # Version: 05/29/2017
  # """)

  t0 = time()

  kappa = 0.2*[1.0 0; 0 1.0]; # conductivity matrix
  magn = 0.06;# heat flux along the boundary
  rin =  1.0;#internal radius

  rex =  2.0; #external radius
  nr = 3; nc = 40;
  Angle = 2*pi;
  thickness =  1.0;
  tolerance = min(rin/nr,  rin/nc/2/pi)/10000;

  fens, fes = Q4annulus(rin, rex, nr, nc, Angle)
  fens, fes = mergenodes(fens,  fes,  tolerance);
  edge_fes = meshboundary(fes);

  # At a single point apply an essential boundary condition (pin down the temperature)
  l1  = selectnode(fens; box=[0.0 0.0 -rex -rex],  inflate = tolerance)
  essential1 = FDataDict("node_list"=>l1, "temperature"=>0.0)

  # The flux boundary condition is applied at two pieces of surface
  # Side 1
  l1 = selectelem(fens, edge_fes, box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
  el1femm = FEMMBase(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)))
  fi = ForceIntensity(FFlt[-magn]);#entering the domain
  flux1 = FDataDict("femm"=>el1femm, "normal_flux"=>-magn) # entering the domain
  # Side 2
  l2=selectelem(fens,edge_fes,box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
  el2femm = FEMMBase(IntegData(subset(edge_fes, l2),  GaussRule(1, 2)))
  flux2 = FDataDict("femm"=>el2femm, "normal_flux"=>+magn) # leaving the domain

  material = MatHeatDiff(kappa)
  femm = FEMMHeatDiff(IntegData(fes,  GaussRule(2, 2)),  material)
  region1 = FDataDict("femm"=>femm)

  # Make model data
  modeldata = FDataDict("fens"=>fens,
  "regions"=>[region1], "essential_bcs"=>[essential1],
  "flux_bcs"=>[flux1, flux2]);

  # Call the solver
  modeldata = FinEtools.AlgoHeatDiffModule.steadystate(modeldata)
  geom=modeldata["geom"]
  Temp=modeldata["temp"]
  # println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")

  # println("Total time elapsed = ",time() - t0,"s")

  # # Postprocessing
  # vtkexportmesh("annulusmod.vtk", fes.conn, [geom.values Temp.values],
  # FinEtools.MeshExportModule.Q8; scalars=[("Temperature", Temp.values)])

  @test abs(minimum(Temp.values)-(-0.50124596))<1.0e-4
  @test abs(maximum(Temp.values)-(+0.50124596))<1.0e-4

  #println("Total time elapsed = ",time() - t0,"s")

  # Postprocessing
  # MeshExportModule.vtkexportmesh ("annulusmod.vtk", fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")
end

end
using .mmmmannulus_Q4_example_algo
mmmmannulus_Q4_example_algo.test()

module mmmmmPoisson_FE_Q4_1
using FinEtools
using Compat.Test
import LinearAlgebra: cholfact
function test()
  # println("""
  #
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit square, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular four-node QUADRILATERALS,
  # in a grid of 1000 x 1000 edges (1M quads, 1M degrees of freedom).
  # Version: 05/29/2017
  # """
  # )
  t0 = time()

  A = 1.0
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = -6.0; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 100;

  # println("Mesh generation")
  fens,fes = Q4block(A, A, N, N)

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))


  # println("Searching nodes  for BC")
  l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
  l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
  l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
  l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
  List = vcat(l1, l2, l3, l4);
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)

  numberdofs!(Temp)

  t1 = time()

  m = MatHeatDiff(thermal_conductivity)
  femm = FEMMHeatDiff(IntegData(fes, GaussRule(2, 2)), m)

  # println("Conductivity")
  K=conductivity(femm, geom, Temp)
  #Profile.print()

  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  fi = ForceIntensity(FFlt, 1, getsource!);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")
  K = cholfact(K)
  # println("Solution of the factorized system")
  U=  K\(F1+F2)
  scattersysvec!(Temp, U[:])


  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  # using .meshExportModule

  # File =  "a.vtk"
  # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")

  Error = 0.0
  for k=1:size(fens.xyz,1)
    Error = Error .+ abs.(Temp.values[k,1] .- tempf(reshape(fens.xyz[k,:], (1,2))))
  end
  # println("Error =$Error")
  @test Error[1]<1.e-5

  true
end
end
using .mmmmmPoisson_FE_Q4_1
mmmmmPoisson_FE_Q4_1.test()

module mmmmmPoisson_FE_example_algo
using FinEtools
using Compat.Test
function test()
    A= 1.0
    thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
    magn = -6.0; #heat source
    truetempf(x) = (1.0 .+ x[1].^2 .+ 2*x[2].^2);#the exact distribution of temperature
    N=20;

    # println("""
    #
    # Heat conduction example described by Amuthan A. Ramabathiran
    # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    # Unit square, with known temperature distribution along the boundary,
    # and uniform heat generation rate inside.  Mesh of regular TRIANGLES,
    # in a grid of $N x $N edges.
    # This version uses the FinEtools algorithm module.
    # """
    # )
    t0 = time()

    fens,fes =T3block(A, A, N, N)


    # Define boundary conditions
    l1 =selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
    l2 =selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
    l3 =selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
    l4 =selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)

    essential1 = FDataDict("node_list"=>vcat(l1, l2, l3, l4),
        "temperature"=>truetempf);
    material = MatHeatDiff(thermal_conductivity)
    femm = FEMMHeatDiff(IntegData(fes, TriRule(1)), material)
    region1 = FDataDict("femm"=>femm, "Q"=>magn)
    # Make model data
    modeldata= FDataDict("fens"=> fens,
        "regions"=>[region1],
        "essential_bcs"=>[essential1]);


    # Call the solver
    modeldata = FinEtools.AlgoHeatDiffModule.steadystate(modeldata)

    # println("Total time elapsed = ",time() - t0,"s")

    geom=modeldata["geom"]
    Temp=modeldata["temp"]
    femm=modeldata["regions"][1]["femm"]
    function errfh(loc,val)
        exact = truetempf(loc)
        return ((exact .- val)^2)[1]
    end

    femm.integdata.integration_rule = TriRule(6)
    E = integratefieldfunction(femm, geom, Temp, errfh, 0.0, m=3)
    # println("Error=$E")

    @test E<00.0025

end
end
using .mmmmmPoisson_FE_example_algo
mmmmmPoisson_FE_example_algo.test()

module mmmmmPoissonRm2
using FinEtools
using Compat.Test
import LinearAlgebra: cholfact
function test()
  # println("""
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit square, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
  # in a grid of 1000 x 1000 edges (2M triangles, 1M degrees of freedom).
  # The material response is defined in a local coordinate system.
  # Version: 05/29/2017
  # """
  # )
  t0 = time()

  A = 1.0 # dimension of the domain (length of the side of the square)
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
  Q = -6.0; # internal heat generation rate
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = Q; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 100;# number of subdivisions along the sides of the square domain
  Rm=[-0.9917568452513019 -0.12813414805267656
  -0.12813414805267656 0.9917568452513019]
  Rm=[-0.8020689950104449 -0.5972313850116512
  -0.5972313850116512 0.8020689950104447]

  # println("Mesh generation")
  fens,fes =T3block(A, A, N, N)

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  # println("Searching nodes  for BC")
  l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
  l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
  l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
  l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
  List = vcat(l1, l2, l3, l4)
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)
  numberdofs!(Temp)

  t1 = time()

  material = MatHeatDiff(thermal_conductivity)

  femm = FEMMHeatDiff(IntegData(fes, TriRule(1)), CSys(Rm), material)


  # println("Conductivity")
  K = conductivity(femm, geom, Temp)
  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  fi = ForceIntensity(FFlt[Q]);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")
  K = cholfact(K)
  # println("Solution of the factorized system")
  U = K\(F1+F2)
  scattersysvec!(Temp,U[:])

  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  Error= 0.0
  for k=1:size(fens.xyz,1)
    Error = Error.+abs.(Temp.values[k,1].-tempf(reshape(fens.xyz[k,:], (1,2))))
  end
  # println("Error =$Error")
  @test Error[1]<1.e-4

end
end
using .mmmmmPoissonRm2
mmmmmPoissonRm2.test()

module mmmmmmmmmNAFEMSm
using FinEtools
using Compat.Test
import LinearAlgebra: norm
function test()
  ## Two-dimensional heat transfer with convection: convergence study
  #

  ## Description
  #
  # Consider a plate of uniform thickness, measuring 0.6 m by 1.0 m. On one
  # short edge the temperature is fixed at 100 °C, and on one long edge the
  # plate is perfectly insulated so that the heat flux is zero through that
  # edge. The other two edges are losing heat via convection to an ambient
  # temperature of 0 °C. The thermal conductivity of the plate is 52.0 W/(m
  # .°K), and the convective heat transfer coefficient is 750 W/(m^2.°K).
  # There is no internal generation of heat. Calculate the temperature 0.2 m
  # along the un-insulated long side, measured from the intersection with the
  # fixed temperature side. The reference result is 18.25 °C.

  ##
  # The reference temperature at the point A  is 18.25 °C according to the
  # NAFEMS publication ( hich cites the book Carslaw, H.S. and J.C. Jaeger,
  # Conduction of Heat in Solids. 1959: Oxford University Press).

  ##
  # The present  tutorial will investigate the reference temperature  and it
  # will attempt to  estimate the  limit value more precisely using a
  # sequence of meshes and Richardson's extrapolation.

  ## Solution
  #

  # println("""
  # NAFEMS benchmark.
  # Two-dimensional heat transfer with convection: convergence study.
  # Solution with quadratic triangles.
  # Version: 05/29/2017
  # """
  # )

  kappa = [52. 0; 0 52.]*phun("W/(M*K)"); # conductivity matrix
  h = 750*phun("W/(M^2*K)");# surface heat transfer coefficient
  Width = 0.6*phun("M");# Geometrical dimensions
  Height = 1.0*phun("M");
  HeightA = 0.2*phun("M");
  Thickness = 0.1*phun("M");
  tolerance  = Width/1000;

  m = MatHeatDiff(kappa)

  modeldata = nothing
  resultsTempA = FFlt[]
  for nref = 1:5
    t0 = time()

    # The mesh is created from two triangles to begin with
    fens,fes = T3blockx([0.0, Width], [0.0, HeightA])
    fens2,fes2 = T3blockx([0.0, Width], [HeightA, Height])
    fens,newfes1,fes2 = mergemeshes(fens, fes, fens2, fes2, tolerance)
    fes = cat(newfes1,fes2)
    # Refine the mesh desired number of times
    for ref = 1:nref
      fens,fes = T3refine(fens,fes);
    end
    fens, fes = T3toT6(fens,fes);
    bfes = meshboundary(fes)

    # Define boundary conditions

    ##
    # The prescribed temperature is applied along edge 1 (the bottom
    # edge in Figure 1)..

    l1 = selectnode(fens; box=[0. Width 0. 0.], inflate=tolerance)
    essential1 = FDataDict("node_list"=>l1, "temperature"=> 100.);

    ##
    # The convection boundary condition is applied along the edges
    # 2,3,4. The elements along the boundary are quadratic line
    # elements L3. The order-four Gauss quadrature is sufficiently
    # accurate.
    l2 = selectelem(fens, bfes; box=[Width Width  0.0 Height], inflate =tolerance)
    l3 = selectelem(fens, bfes; box=[0.0 Width Height Height], inflate =tolerance)
    cfemm = FEMMHeatDiffSurf(IntegData(subset(bfes,vcat(l2,l3)),
      GaussRule(1, 3), Thickness), h)
    convection1 = FDataDict("femm"=>cfemm, "ambient_temperature"=>0.);

    # The interior
    femm = FEMMHeatDiff(IntegData(fes, TriRule(3), Thickness), m)
    region1 = FDataDict("femm"=>femm)

    # Make the model data
    modeldata = FDataDict("fens"=> fens,
    "regions"=>[region1],
    "essential_bcs"=>[essential1],
    "convection_bcs"=>[convection1]);

    # Call the solver
    modeldata = FinEtools.AlgoHeatDiffModule.steadystate(modeldata)

    # println("Total time elapsed = ",time() - t0,"s")

    l4 = selectnode(fens; box=[Width Width HeightA HeightA], inflate =tolerance)

    geom = modeldata["geom"]
    Temp = modeldata["temp"]

    ##
    # Collect the temperature  at the point A  [coordinates
    # (Width,HeightA)].
    push!(resultsTempA, Temp.values[l4][1]);

  end

  ##
  # These are the computed results for the temperature at point A:
  # println("$( resultsTempA  )")

  # Postprocessing
  geom = modeldata["geom"]
  Temp = modeldata["temp"]
  regions = modeldata["regions"]
  vtkexportmesh("T4NAFEMS--T6.vtk", connasarray(regions[1]["femm"].integdata.fes),
  [geom.values Temp.values/100], FinEtools.MeshExportModule.T6;
  scalars=[("Temperature", Temp.values)])
  try  rm("T4NAFEMS--T6.vtk"); catch end
  vtkexportmesh("T4NAFEMS--T6--base.vtk", connasarray(regions[1]["femm"].integdata.fes),
  [geom.values 0.0*Temp.values/100], FinEtools.MeshExportModule.T6)
  try rm("T4NAFEMS--T6--base.vtk"); catch end
  # ##
  # # Richardson extrapolation is used to estimate the true solution from the
  # # results for the finest three meshes.
  #    [xestim, beta] = richextrapol(results(end-2:end),mesh_sizes(end-2:end));
  #     disp(['Estimated true solution for temperature at A: ' num2str(xestim) ' degrees'])

  # ##
  # # Plot the estimated true error.
  #    figure
  #     loglog(mesh_sizes,abs(results-xestim)/xestim,'bo-','linewidth',3)
  #     grid on
  #      xlabel('log(mesh size)')
  #     ylabel('log(|estimated temperature error|)')
  #     set_graphics_defaults

  # ##
  # # The estimated true error has  a slope of approximately 4 on the log-log
  # scale.
  # ##
  # # Plot the absolute values of the approximate error (differences  of
  # # successive solutions).
  #     figure
  #     loglog(mesh_sizes(2:end),abs(diff(results)),'bo-','linewidth',3)
  #     Thanksgrid on
  #     xlabel('log(mesh size)')
  #     ylabel('log(|approximate temperature error|)')
  #     set_graphics_defaults


  ## Discussion
  #
  ##
  # The last segment  of the approximate error curve is close to the slope of
  # the estimated true error. Nevertheless, it would have been more
  # reassuring if the  three successive approximate errors  were located more
  # closely on a straight line.

  ##
  # The use of uniform mesh-size meshes is sub optimal: it would be more
  # efficient to use graded meshes. The tutorial pub_T4NAFEMS_conv_graded
  # addresses use of graded meshes  in convergence studies.


  @test (norm(resultsTempA-[17.9028, 18.3323, 18.2965, 18.2619, 18.255]))<1.0e-3

end
end
using .mmmmmmmmmNAFEMSm
mmmmmmmmmNAFEMSm.test()

module mmmmmmconvergence
using FinEtools
using Compat.Test
import LinearAlgebra: norm
function test()
  ## Two-dimensional heat transfer with convection: convergence study
  #

  ## Description
  #
  # Consider a plate of uniform thickness, measuring 0.6 m by 1.0 m. On one
  # short edge the temperature is fixed at 100 °C, and on one long edge the
  # plate is perfectly insulated so that the heat flux is zero through that
  # edge. The other two edges are losing heat via convection to an ambient
  # temperature of 0 °C. The thermal conductivity of the plate is 52.0 W/(m
  # .°K), and the convective heat transfer coefficient is 750 W/(m^2.°K).
  # There is no internal generation of heat. Calculate the temperature 0.2 m
  # along the un-insulated long side, measured from the intersection with the
  # fixed temperature side. The reference result is 18.25 °C.

  ##
  # The reference temperature at the point A  is 18.25 °C according to the
  # NAFEMS publication ( hich cites the book Carslaw, H.S. and J.C. Jaeger,
  # Conduction of Heat in Solids. 1959: Oxford University Press).

  ##
  # The present  tutorial will investigate the reference temperature  and it
  # will attempt to  estimate the  limit value more precisely using a
  # sequence of meshes and Richardson's extrapolation.

  ## Solution
  #

  # println("""
  # NAFEMS benchmark.
  # Two-dimensional heat transfer with convection: convergence study.
  # Solution with linear triangles.
  # Version: 05/29/2017
  # """
  # )

  kappa = [52. 0; 0 52.]*phun("W/(M*K)"); # conductivity matrix
  h = 750*phun("W/(M^2*K)");# surface heat transfer coefficient
  Width = 0.6*phun("M");# Geometrical dimensions
  Height = 1.0*phun("M");
  HeightA = 0.2*phun("M");
  Thickness = 0.1*phun("M");
  tolerance  = Width/1000;

  m = MatHeatDiff(kappa)

  modeldata = nothing
  resultsTempA = FFlt[]
  for nref = 1:5
    t0 = time()

    # The mesh is created from two triangles to begin with
    fens,fes = T3blockx([0.0, Width], [0.0, HeightA])
    fens2,fes2 = T3blockx([0.0, Width], [HeightA, Height])
    fens,newfes1,fes2 = mergemeshes(fens, fes, fens2, fes2, tolerance)
    fes = cat(newfes1,fes2)
    # Refine the mesh desired number of times
    for ref = 1:nref
      fens,fes = T3refine(fens,fes);
    end
    bfes = meshboundary(fes)

    # Define boundary conditions

    ##
    # The prescribed temperature is applied along edge 1 (the bottom
    # edge in Figure 1)..

    l1 = selectnode(fens; box=[0. Width 0. 0.], inflate=tolerance)
    essential1 = FDataDict("node_list"=>l1, "temperature"=> 100.);

    ##
    # The convection boundary condition is applied along the edges
    # 2,3,4. The elements along the boundary are quadratic line
    # elements L3. The order-four Gauss quadrature is sufficiently
    # accurate.
    l2 = selectelem(fens, bfes; box=[Width Width  0.0 Height], inflate =tolerance)
    l3 = selectelem(fens, bfes; box=[0.0 Width Height Height], inflate =tolerance)
    cfemm = FEMMHeatDiffSurf(IntegData(subset(bfes,vcat(l2,l3)),
      GaussRule(1, 3), Thickness), h)
    convection1 = FDataDict("femm"=>cfemm, "ambient_temperature"=>0.);

    # The interior
    femm = FEMMHeatDiff(IntegData(fes, TriRule(3), Thickness), m)
    region1 = FDataDict("femm"=>femm)

    # Make the model data
    modeldata = FDataDict("fens"=> fens,
    "regions"=>[region1],
    "essential_bcs"=>[essential1],
    "convection_bcs"=>[convection1]);

    # Call the solver
    modeldata = FinEtools.AlgoHeatDiffModule.steadystate(modeldata)

    # println("Total time elapsed = ",time() - t0,"s")

    l4 = selectnode(fens; box=[Width Width HeightA HeightA], inflate =tolerance)

    geom = modeldata["geom"]
    Temp = modeldata["temp"]

    ##
    # Collect the temperature  at the point A  [coordinates
    # (Width,HeightA)].
    push!(resultsTempA, Temp.values[l4][1]);

  end

  ##
  # These are the computed results for the temperature at point A:
  # println("$( resultsTempA  )")

  # Postprocessing
  geom = modeldata["geom"]
  Temp = modeldata["temp"]
  regions = modeldata["regions"]
  vtkexportmesh("T4NAFEMS--T3.vtk", connasarray(regions[1]["femm"].integdata.fes),
  [geom.values Temp.values/100], FinEtools.MeshExportModule.T3;
  scalars=[("Temperature", Temp.values)])
  rm("T4NAFEMS--T3.vtk")
  vtkexportmesh("T4NAFEMS--T3--base.vtk", connasarray(regions[1]["femm"].integdata.fes),
  [geom.values 0.0*Temp.values/100], FinEtools.MeshExportModule.T3)
  rm("T4NAFEMS--T3--base.vtk")
  # ##
  # # Richardson extrapolation is used to estimate the true solution from the
  # # results for the finest three meshes.
  #    [xestim, beta] = richextrapol(results(end-2:end),mesh_sizes(end-2:end));
  #     disp(['Estimated true solution for temperature at A: ' num2str(xestim) ' degrees'])

  # ##
  # # Plot the estimated true error.
  #    figure
  #     loglog(mesh_sizes,abs(results-xestim)/xestim,'bo-','linewidth',3)
  #     grid on
  #      xlabel('log(mesh size)')
  #     ylabel('log(|estimated temperature error|)')
  #     set_graphics_defaults

  # ##
  # # The estimated true error has  a slope of approximately 4 on the log-log
  # scale.
  # ##
  # # Plot the absolute values of the approximate error (differences  of
  # # successive solutions).
  #     figure
  #     loglog(mesh_sizes(2:end),abs(diff(results)),'bo-','linewidth',3)
  #     Thanksgrid on
  #     xlabel('log(mesh size)')
  #     ylabel('log(|approximate temperature error|)')
  #     set_graphics_defaults


  ## Discussion
  #
  ##
  # The last segment  of the approximate error curve is close to the slope of
  # the estimated true error. Nevertheless, it would have been more
  # reassuring if the  three successive approximate errors  were located more
  # closely on a straight line.

  ##
  # The use of uniform mesh-size meshes is sub optimal: it would be more
  # efficient to use graded meshes. The tutorial pub_T4NAFEMS_conv_graded
  # addresses use of graded meshes  in convergence studies.

  @test (norm(resultsTempA- [22.7872, 19.1813, 18.516, 18.3816, 18.3064])       )<1.0e-3

end
end
using .mmmmmmconvergence
mmmmmmconvergence.test()

module mmmPoissonmmmt4mmmm
using FinEtools
using Compat.Test
function test()

  t0 = time()

  A = 1.0 # dimension of the domain (length of the side of the square)
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]; # conductivity matrix
  Q = -6.0; # internal heat generation rate
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = Q; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 30;# number of subdivisions along the sides of the square domain

  # println("Mesh generation")
  fens,fes = T4block(A, A, A, N, N, N)

  # println("""
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit cube, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular quadratic TETRAHEDRA,
  # in a grid of $(N) x $(N) x $(N) edges ($(count(fens)) degrees of freedom).
  # Version: 07/03/2017
  # """
  # )

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  # println("Searching nodes  for BC")
  Tolerance = 1.0/N/100.0
  l1 = selectnode(fens; box=[0. 0. 0. A 0. A], inflate = Tolerance)
  l2 = selectnode(fens; box=[A A 0. A 0. A], inflate = Tolerance)
  l3 = selectnode(fens; box=[0. A 0. 0. 0. A], inflate = Tolerance)
  l4 = selectnode(fens; box=[0. A A A 0. A], inflate = Tolerance)
  l5 = selectnode(fens; box=[0. A 0. A 0. 0.], inflate = Tolerance)
  l6 = selectnode(fens; box=[0. A 0. A A A], inflate = Tolerance)
  List = vcat(l1, l2, l3, l4, l5, l6)
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)
  numberdofs!(Temp)

  # println( "Number of free degrees of freedom: $(Temp.nfreedofs)")
  t1 = time()

  material = MatHeatDiff(thermal_conductivity)

  femm = FEMMHeatDiff(IntegData(fes, TetRule(1), 100.), material)


  # println("Conductivity")
  K = conductivity(femm, geom, Temp)
  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
  fi = ForceIntensity(FFlt[Q]);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")


  # println("Solution of the factorized system")
  U = K\(F1+F2)
  scattersysvec!(Temp,U[:])
  #
  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  Error= 0.0
  for k=1:size(fens.xyz,1)
    Error = Error.+abs.(Temp.values[k,1].-tempf(reshape(fens.xyz[k,:], (1,3))))
  end
  # println("Error =$Error")
@test abs(Error[1]) < 1.0e-3
end
end
using .mmmPoissonmmmt4mmmm
mmmPoissonmmmt4mmmm.test()

module mmmPoissonmmmt10mm
using FinEtools
using Compat.Test
function test()

  t0 = time()

  A = 1.0 # dimension of the domain (length of the side of the square)
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]; # conductivity matrix
  Q = -6.0; # internal heat generation rate
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = Q; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 13;# number of subdivisions along the sides of the square domain

  # println("Mesh generation")
  fens,fes = T10block(A, A, A, N, N, N)

  # println("""
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit cube, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular quadratic TETRAHEDRA,
  # in a grid of $(N) x $(N) x $(N) edges ($(count(fens)) degrees of freedom).
  # Version: 07/03/2017
  # """
  # )

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  # println("Searching nodes  for BC")
  Tolerance = 1.0/N/100.0
  l1 = selectnode(fens; box=[0. 0. 0. A 0. A], inflate = Tolerance)
  l2 = selectnode(fens; box=[A A 0. A 0. A], inflate = Tolerance)
  l3 = selectnode(fens; box=[0. A 0. 0. 0. A], inflate = Tolerance)
  l4 = selectnode(fens; box=[0. A A A 0. A], inflate = Tolerance)
  l5 = selectnode(fens; box=[0. A 0. A 0. 0.], inflate = Tolerance)
  l6 = selectnode(fens; box=[0. A 0. A A A], inflate = Tolerance)
  List = vcat(l1, l2, l3, l4, l5, l6)
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)
  numberdofs!(Temp)

  # println( "Number of free degrees of freedom: $(Temp.nfreedofs)")
  t1 = time()

  material = MatHeatDiff(thermal_conductivity)

  femm = FEMMHeatDiff(IntegData(fes, TetRule(4), 100.), material)


  # println("Conductivity")
  K = conductivity(femm, geom, Temp)
  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
  fi = ForceIntensity(FFlt[Q]);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")


  # println("Solution of the factorized system")
  U = K\(F1+F2)
  scattersysvec!(Temp,U[:])
  #
  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  Error= 0.0
  for k=1:size(fens.xyz,1)
    Error = Error .+ abs.(Temp.values[k,1] .- tempf(reshape(fens.xyz[k,:], (1,3))))
  end
  # println("Error =$Error")
@test abs(Error[1]) < 1.0e-3
end
end
using .mmmPoissonmmmt10mm
mmmPoissonmmmt10mm.test()

module mmmPoissonmmh8mm
using FinEtools
using Compat.Test
function test()

  t0 = time()

  A = 1.0 # dimension of the domain (length of the side of the square)
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]; # conductivity matrix
  Q = -6.0; # internal heat generation rate
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = Q; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 13;# number of subdivisions along the sides of the square domain

  # println("Mesh generation")
  fens,fes = H8block(A, A, A, N, N, N)

  # println("""
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit cube, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular quadratic TETRAHEDRA,
  # in a grid of $(N) x $(N) x $(N) edges ($(count(fens)) degrees of freedom).
  # Version: 07/03/2017
  # """
  # )

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  # println("Searching nodes  for BC")
  Tolerance = 1.0/N/100.0
  l1 = selectnode(fens; box=[0. 0. 0. A 0. A], inflate = Tolerance)
  l2 = selectnode(fens; box=[A A 0. A 0. A], inflate = Tolerance)
  l3 = selectnode(fens; box=[0. A 0. 0. 0. A], inflate = Tolerance)
  l4 = selectnode(fens; box=[0. A A A 0. A], inflate = Tolerance)
  l5 = selectnode(fens; box=[0. A 0. A 0. 0.], inflate = Tolerance)
  l6 = selectnode(fens; box=[0. A 0. A A A], inflate = Tolerance)
  List = vcat(l1, l2, l3, l4, l5, l6)
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)
  numberdofs!(Temp)

  # println( "Number of free degrees of freedom: $(Temp.nfreedofs)")
  t1 = time()

  material = MatHeatDiff(thermal_conductivity)

  femm = FEMMHeatDiff(IntegData(fes, GaussRule(3,2), 100.), material)


  # println("Conductivity")
  K = conductivity(femm, geom, Temp)
  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
  fi = ForceIntensity(FFlt[Q]);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")


  # println("Solution of the factorized system")
  U = K\(F1+F2)
  scattersysvec!(Temp,U[:])
  #
  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  Error= 0.0
  for k=1:size(fens.xyz,1)
    Error = Error .+ abs.(Temp.values[k,1] .- tempf(reshape(fens.xyz[k,:], (1,3))))
  end
  # println("Error =$Error")
@test abs(Error[1]) < 1.0e-3
end
end
using .mmmPoissonmmh8mm
mmmPoissonmmh8mm.test()

module mmmmmactuatormmh8m
using FinEtools
using Compat.Test
function test()
  # MEMS actuator.   Thermal analysis.
  x0 =  0.0*phun("micro*m");
  x1 = x0+5.0*phun("micro*m");
  x2 = x1+10.0*phun("micro*m");
  x3 = x2+10.0*phun("micro*m");
  x4 = x3+10.0*phun("micro*m");
  y0 = 0.0*phun("micro*m");
  y4 = 250.0*phun("micro*m");
  y3 = y4-10.0*phun("micro*m");
  y2 = y3-10.0*phun("micro*m");
  y1 = y2-10.0*phun("micro*m");
  t = 2.0*phun("micro*m");
  h = 0.1*phun("micro*m");
  z0 = 0.0*phun("micro*m");
  z3 = 2*t+h;
  z2 = z3-t;
  z1 = z2-h;
  m1 = 2*2;
  m2 = 2*2;
  m3 = 2*2;
  m4 = 3*2;
  n1 = 20*2;
  n2 = 4*2;
  n3 = 2*2;
  n4 = 2*2;
  n5 = 7*2;
  p1 = 1*2;
  p2 = 1*2;
  p3 = 1*2;
  kappa = 157*[i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]*phun("W/m/K"); # W/m/K, conductivity matrix
  DV = 5*phun("V"); # voltage drop in volt
  l  = 2*(y1+y2)/2+2*(x1+x2)/2; # length of the conductor
  resistivity  =  1.1e-5*phun("Ohm*m"); # Ohm m
  Q = DV^2/resistivity/l^2; # rate of Joule heating, W/m^3
  T_substrate = 293; # substrate temperature in degrees Kelvin

  fens,fes =  H8hexahedron([x1 y0 z0; x2 y1 z1],m2,n1,p1);
  fens1,fes1  =  H8hexahedron([x1 y1 z0;x2 y2 z1],m2,n2,p1);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x0 y1 z0;x1 y2 z1],m1,n2,p1);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x0 y1 z1;x1 y2 z2], m1,n2,p2);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x0 y1 z2;x1 y2 z3],m1,n2,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x0 y2 z2;x1 y3 z3],m1,n3,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x0 y3 z2;x1 y4 z3], m1,n4,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x1 y3 z2;x3 y4 z3],m4,n4,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x3 y3 z2;x4 y4 z3],m3,n4,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x3 y0 z2;x4 y3 z3], m3,n5,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);

  hotmater = MatHeatDiff(kappa)
  coldmater = MatHeatDiff(kappa)
   cl =  selectelem(fens, fes, box=[x0,x2,y0,y2,z0,z1],inflate = t/100);

  hotfemm  =  FEMMHeatDiff(IntegData(subset(fes,cl), GaussRule(3, 2), 0.), hotmater)
  coldfemm  = FEMMHeatDiff(IntegData(subset(fes,setdiff(collect(1:count(fes)), cl)),
    GaussRule(3, 2), 0.), coldmater)
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
  fenids = selectnode(fens, box=[x0,x4,y0,y0,z0,z3],
      inflate=t/1000) ; # fixed temperature on substrate
  setebc!(Temp, fenids, true, 1, T_substrate)
  applyebc!(Temp)
  numberdofs!(Temp)

  K = conductivity(hotfemm, geom, Temp) +
  conductivity(coldfemm, geom, Temp)
  fi = ForceIntensity(FFlt[Q]);
  F = distribloads(hotfemm, geom, Temp, fi, 3) +
  nzebcloadsconductivity(hotfemm, geom, Temp) +
  nzebcloadsconductivity(coldfemm, geom, Temp)

  U = K\F
  scattersysvec!(Temp,U[:])

  # using Plots
  # plotly()
  nList = selectnode(fens, box=[x1,x1,y0,y1,z1,z1], inflate=t/100)
  y_i = geom.values[nList, 2]
  T_i = Temp.values[nList, 1]
  ix = sortperm(y_i)
  # plot(y_i[ix], T_i[ix], color=:red, label= "hot leg")

  nList = selectnode(fens, box=[x3,x3,y0,y3,z2,z2], inflate=t/100)
  y_o = geom.values[nList, 2]
  T_o = Temp.values[nList, 1]
  ix = sortperm(y_o)
  # plot!(y_o[ix], T_o[ix], color=:blue, label= "cold leg")

# show(T_i)
# println("maximum(T_i) = $(maximum(T_i))")
@test abs(maximum(T_i)-1380.5883006341187) < 1.0e-3
end
end
using .mmmmmactuatormmh8m
mmmmmactuatormmh8m.test()

module mmmmmactuatormmmmmmmm
using FinEtools
using Compat.Test
function test()
  # MEMS actuator.   Thermal analysis.
  x0 =  0.0*phun("micro*m");
  x1 = x0+5.0*phun("micro*m");
  x2 = x1+10.0*phun("micro*m");
  x3 = x2+10.0*phun("micro*m");
  x4 = x3+10.0*phun("micro*m");
  y0 = 0.0*phun("micro*m");
  y4 = 250.0*phun("micro*m");
  y3 = y4-10.0*phun("micro*m");
  y2 = y3-10.0*phun("micro*m");
  y1 = y2-10.0*phun("micro*m");
  t = 2.0*phun("micro*m");
  h = 0.1*phun("micro*m");
  z0 = 0.0*phun("micro*m");
  z3 = 2*t+h;
  z2 = z3-t;
  z1 = z2-h;
  m1 = 2*2;
  m2 = 2*2;
  m3 = 2*2;
  m4 = 3*2;
  n1 = 20*2;
  n2 = 4*2;
  n3 = 2*2;
  n4 = 2*2;
  n5 = 7*2;
  p1 = 1*2;
  p2 = 1*2;
  p3 = 1*2;
  kappa = 157*[i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]*phun("W/m/K"); # W/m/K, conductivity matrix
  DV = 5*phun("V"); # voltage drop in volt
  l  = 2*(y1+y2)/2+2*(x1+x2)/2; # length of the conductor
  resistivity  =  1.1e-5*phun("Ohm*m"); # Ohm m
  Q = DV^2/resistivity/l^2; # rate of Joule heating, W/m^3
  T_substrate = 293; # substrate temperature in degrees Kelvin

  fens,fes =  H8hexahedron([x1 y0 z0; x2 y1 z1],m2,n1,p1);
  fens1,fes1  =  H8hexahedron([x1 y1 z0;x2 y2 z1],m2,n2,p1);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x0 y1 z0;x1 y2 z1],m1,n2,p1);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x0 y1 z1;x1 y2 z2], m1,n2,p2);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x0 y1 z2;x1 y2 z3],m1,n2,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x0 y2 z2;x1 y3 z3],m1,n3,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x0 y3 z2;x1 y4 z3], m1,n4,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x1 y3 z2;x3 y4 z3],m4,n4,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x3 y3 z2;x4 y4 z3],m3,n4,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);
  fens1,fes1  =  H8hexahedron([x3 y0 z2;x4 y3 z3], m3,n5,p3);
  fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
  fes =  cat(fes1,fes2);

  fens,fes  =  H8toH20(fens,fes);

  hotmater = MatHeatDiff(kappa)
  coldmater = MatHeatDiff(kappa)
   cl =  selectelem(fens, fes, box=[x0,x2,y0,y2,z0,z1],inflate = t/100);

  hotfemm  =  FEMMHeatDiff(IntegData(subset(fes,cl), GaussRule(3, 3), 0.), hotmater)
  coldfemm  = FEMMHeatDiff(IntegData(subset(fes,setdiff(collect(1:count(fes)), cl)),
    GaussRule(3, 3), 0.), coldmater)
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
  fenids = selectnode(fens, box=[x0,x4,y0,y0,z0,z3],
      inflate=t/1000) ; # fixed temperature on substrate
  setebc!(Temp, fenids, true, 1, T_substrate)
  applyebc!(Temp)
  numberdofs!(Temp)

  K = conductivity(hotfemm, geom, Temp) +
  conductivity(coldfemm, geom, Temp)
  fi = ForceIntensity(FFlt[Q]);
  F = distribloads(hotfemm, geom, Temp, fi, 3) +
  nzebcloadsconductivity(hotfemm, geom, Temp) +
  nzebcloadsconductivity(coldfemm, geom, Temp)

  U = K\F
  scattersysvec!(Temp,U[:])

  # using Plots
  # plotly()
  nList = selectnode(fens, box=[x1,x1,y0,y1,z1,z1], inflate=t/100)
  y_i = geom.values[nList, 2]
  T_i = Temp.values[nList, 1]
  ix = sortperm(y_i)
  # plot(y_i[ix], T_i[ix], color=:red, label= "hot leg")

  nList = selectnode(fens, box=[x3,x3,y0,y3,z2,z2], inflate=t/100)
  y_o = geom.values[nList, 2]
  T_o = Temp.values[nList, 1]
  ix = sortperm(y_o)
  # plot!(y_o[ix], T_o[ix], color=:blue, label= "cold leg")

# show(T_i)
# println("maximum(T_i) = $(maximum(T_i))")
@test abs(maximum(T_i)-1382.089817533287) < 1.0e-3
end
end
using .mmmmmactuatormmmmmmmm
mmmmmactuatormmmmmmmm.test()

module mmmmPoiss_07082017
using FinEtools
using Compat.Test
function test()

  # println("""
  #
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit square, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
  # in a grid of 1000 x 1000 edges (2M triangles, 1M degrees of freedom).
  # Version: 05/29/2017
  # """
  # )
  t0 = time()

  A = 1.0 # dimension of the domain (length of the side of the square)
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
  Q = -6.0; # internal heat generation rate
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = Q; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 100;# number of subdivisions along the sides of the square domain


  # println("Mesh generation")
  fens,fes =T3block(A, A, N, N)

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  # println("Searching nodes  for BC")
  l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
  l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
  l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
  l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
  List = vcat(l1, l2, l3, l4)
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)
  numberdofs!(Temp)

  t1 = time()

  material = MatHeatDiff(thermal_conductivity)

  femm = FEMMHeatDiff(IntegData(fes, TriRule(1)), CSys(2, 2), material)


  # println("Conductivity")
  K = conductivity(femm, geom, Temp)
  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  fi = ForceIntensity(FFlt, 1, getsource!);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")
  # K = cholfact(K)
  # println("Solution of the factorized system")
  U = K\(F1+F2)
  scattersysvec!(Temp,U[:])

  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  Error= 0.0
  for k=1:size(fens.xyz,1)
    Error = Error .+ abs.(Temp.values[k,1] .- tempf(reshape(fens.xyz[k,:], (1,2))))
  end
  # println("Error =$Error")


  # File =  "a.vtk"
  # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

  @test Error[1]<1.e-5

  true
end
end
using .mmmmPoiss_07082017
mmmmPoiss_07082017.test()

module mmmmmAnnularQ8mmmmm
using FinEtools
using Compat.Test
import LinearAlgebra: norm, cholfact
function test()

  # println("""
  # Annular region,  ingoing and outgoing flux. Minimum/maximum temperature ~(+/-)0.501.
  # Mesh of quadratic serendipity quadrilaterals.
  # Version: 05/29/2017
  # """)

  t0 = time()

  kappa = 0.2*[1.0 0; 0 1.0]; # conductivity matrix
  magn = 0.06;# heat flux along the boundary
  rin =  1.0;#internal radius

  rex =  2.0; #external radius
  nr = 5; nc = 40;
  Angle = 2*pi;
  thickness =  1.0;
  tolerance = min(rin/nr,  rin/nc/2/pi)/10000;


  fens, fes  =  Q8annulus(rin, rex, nr, nc, Angle)
  fens, fes  =  mergenodes(fens,  fes,  tolerance);
  edge_fes  =  meshboundary(fes);

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz, 1), 1))


  l1  = selectnode(fens; box=[0.0 0.0 -rex -rex],  inflate = tolerance)
  setebc!(Temp, l1, 1; val=zero(FFlt))
  applyebc!(Temp)

  numberdofs!(Temp)


  material = MatHeatDiff(kappa)
  femm = FEMMHeatDiff(IntegData(fes,  GaussRule(2, 3)),  material)

  K = conductivity(femm,  geom,  Temp)

  l1 = selectelem(fens, edge_fes, box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
  el1femm = FEMMBase(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)))
  fi = ForceIntensity(FFlt[-magn]);#entering the domain
  F1 = (-1.0)* distribloads(el1femm,  geom,  Temp,  fi,  2);

  l1 = selectelem(fens, edge_fes, box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
  el1femm =  FEMMBase(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)))
  fi = ForceIntensity(FFlt[+magn]);#leaving the domain
  F2 = (-1.0)* distribloads(el1femm,  geom,  Temp,  fi,  2);

  F3 = nzebcloadsconductivity(femm,  geom,  Temp);


  K = cholfact(K)
  U = K\(F1+F2+F3)
  scattersysvec!(Temp, U[:])

  # println("Total time elapsed = ", time() - t0, "s")

  File =  "annulusq8.vtk"
  vtkexportmesh(File,  connasarray(fes),  [geom.values Temp.values],
  FinEtools.MeshExportModule.Q8; scalars=[("Temperature", Temp.values)])
  try rm(File); catch end
  # println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")
  @test norm([minimum(Temp.values), maximum(Temp.values)]-[-0.5010001850658392, 0.5010001850658563]) < 1.0e-5
  true

end
end
using .mmmmmAnnularQ8mmmmm
mmmmmAnnularQ8mmmmm.test()

module mmpnpConcreteColumnsymb_conv
using FinEtools
using Compat.Test
import LinearAlgebra: norm, cholfact
function test()
  a=2.5; dy=a/2*sin(15/180*pi); dx=a/2*cos(15/180*pi); Q=4.5; k=1.8; Dz=1.0;
  h= 5.;

  m = MatHeatDiff([k 0; 0 k])

  modeldata = nothing

  fens = FENodeSet([0 0; dx -dy; dx dy; 2*dx -2*dy; 2*dx 2*dy])
  fes = FESetT3([1 2 3; 2 4 5; 2 5 3]);

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  applyebc!(Temp)
  numberdofs!(Temp)
  Temp.dofnums[1] = 3
  Temp.dofnums[2] = 1
  Temp.dofnums[3] = 2
  Temp.dofnums[4] = 5
  Temp.dofnums[5] = 4

  bfes = FESetL2([4 5]);

  cfemm = FEMMHeatDiffSurf(IntegData(bfes, GaussRule(1, 2), Dz), h)
  femm = FEMMHeatDiff(IntegData(fes, TriRule(1), Dz), m)
  fi = ForceIntensity(FFlt[Q]);
  F1 = distribloads(femm, geom, Temp, fi, 3);
  K = conductivity(femm, geom, Temp)
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  H = surfacetransfer(cfemm, geom, Temp);
  # show(H)
  # F3 = surfacetransferloads(femm, geom, temp, amb);
  Factor = cholfact(K+H)
  U = Factor\(F1+F2)
  scattersysvec!(Temp, U[:])
  # display(Temp)
  @test norm([5.1362
              3.98612
              3.85657
              1.00778
              1.16556] - Temp.values) < 1.0e-3
end
end
using .mmpnpConcreteColumnsymb_conv
mmpnpConcreteColumnsymb_conv.test()

module mmmmPoiss_08082017
using FinEtools
using FinEtools.MeshExportModule
using Compat.Test
import LinearAlgebra: cholfact
function test()

  # println("""
  #
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit square, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
  # in a grid of 1000 x 1000 edges (2M triangles, 1M degrees of freedom).
  # Version: 05/29/2017
  # """
  # )
  t0 = time()

  A = 1.0 # dimension of the domain (length of the side of the square)
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
  Q = -6.0; # internal heat generation rate
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = Q; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 10;# number of subdivisions along the sides of the square domain


  # println("Mesh generation")
  fens,fes =T3block(A, A, N, N)
  fens,fes = T3toT6(fens,fes)

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  # println("Searching nodes  for BC")
  l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
  l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
  l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
  l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
  List = vcat(l1, l2, l3, l4)
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)
  numberdofs!(Temp)

  t1 = time()

  material = MatHeatDiff(thermal_conductivity)

  femm = FEMMHeatDiff(IntegData(fes, TriRule(13), 100.), material)


  # println("Conductivity")
  K = conductivity(femm, geom, Temp)
  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
  fi = ForceIntensity(FFlt[Q]);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")
  K = cholfact(K)
  # println("Solution of the factorized system")
  U = K\(F1+F2)
  scattersysvec!(Temp,U[:])

  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  Error= 0.0
  for k=1:size(fens.xyz,1)
    Error = Error .+ abs.(Temp.values[k,1] .- tempf(reshape(fens.xyz[k,:], (1,2))))
  end
  # println("Error =$Error")


  # File =  "a.vtk"
  # MeshExportModule.vtkexportmesh(File, fes.conn, hcat(geom.values, Temp.values), MeshExportModule.T6;
  # scalars=("Temperature", Temp.values))

  @test Error[1]<1.e-5

  true
end
end
using .mmmmPoiss_08082017
mmmmPoiss_08082017.test()


module mmmmmPoisson_FE_Q8_2
using FinEtools
using FinEtools.MeshExportModule
using Compat.Test
import LinearAlgebra: cholfact
function test()
  # println("""
  #
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit square, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular four-node QUADRILATERALS,
  # in a grid of 1000 x 1000 edges (1M quads, 1M degrees of freedom).
  # Version: 05/29/2017
  # """
  # )
  t0 = time()

  A = 1.0
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = -6.0; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 10;

  # println("Mesh generation")
  fens,fes = Q8block(A, A, N, N)

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))


  # println("Searching nodes  for BC")
  l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
  l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
  l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
  l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
  List = vcat(l1, l2, l3, l4);
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)

  numberdofs!(Temp)

  t1 = time()

  m = MatHeatDiff(thermal_conductivity)
  femm = FEMMHeatDiff(IntegData(fes, GaussRule(2, 4)), m)

  # println("Conductivity")
  K=conductivity(femm, geom, Temp)
  #Profile.print()

  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  fi = ForceIntensity(FFlt, 1, getsource!);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")
  K = cholfact(K)
  # println("Solution of the factorized system")
  U=  K\(F1+F2)
  scattersysvec!(Temp, U[:])


  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  # using .meshExportModule

  # File =  "a.vtk"
  # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")

  Error = 0.0
  for k=1:size(fens.xyz,1)
    Error = Error .+ abs.(Temp.values[k,1] .- tempf(reshape(fens.xyz[k,:], (1,2))))
  end
  # println("Error =$Error")
  @test Error[1]<1.e-5

    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh(File, fes.conn, hcat(geom.values, Temp.values),
    # MeshExportModule.Q8;
    # scalars=("Temperature", Temp.values))

  true
end
end
using .mmmmmPoisson_FE_Q8_2
mmmmmPoisson_FE_Q8_2.test()


module mmmPoissonmmmt10m2m
using FinEtools
using Compat.Test
function test()

  t0 = time()

  A = 1.0 # dimension of the domain (length of the side of the square)
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]; # conductivity matrix
  Q = -6.0; # internal heat generation rate
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = Q; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 10;# number of subdivisions along the sides of the square domain

  # println("Mesh generation")
  fens,fes = T10block(A, A, A, N, N, N)

  # println("""
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit cube, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular quadratic TETRAHEDRA,
  # in a grid of $(N) x $(N) x $(N) edges ($(count(fens)) degrees of freedom).
  # Version: 07/03/2017
  # """
  # )

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  # println("Searching nodes  for BC")
  Tolerance = 1.0/N/100.0
  l1 = selectnode(fens; box=[0. 0. 0. A 0. A], inflate = Tolerance)
  l2 = selectnode(fens; box=[A A 0. A 0. A], inflate = Tolerance)
  l3 = selectnode(fens; box=[0. A 0. 0. 0. A], inflate = Tolerance)
  l4 = selectnode(fens; box=[0. A A A 0. A], inflate = Tolerance)
  l5 = selectnode(fens; box=[0. A 0. A 0. 0.], inflate = Tolerance)
  l6 = selectnode(fens; box=[0. A 0. A A A], inflate = Tolerance)
  List = vcat(l1, l2, l3, l4, l5, l6)
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)
  numberdofs!(Temp)

  # println( "Number of free degrees of freedom: $(Temp.nfreedofs)")
  t1 = time()

  material = MatHeatDiff(thermal_conductivity)

  femm = FEMMHeatDiff(IntegData(fes, TetRule(5), 100.), material)


  # println("Conductivity")
  K = conductivity(femm, geom, Temp)
  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
  fi = ForceIntensity(FFlt[Q]);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")


  # println("Solution of the factorized system")
  U = K\(F1+F2)
  scattersysvec!(Temp,U[:])
  #
  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  Error= 0.0
  for k=1:size(fens.xyz,1)
    Error = Error .+ abs.(Temp.values[k,1] .- tempf(reshape(fens.xyz[k,:], (1,3))))
  end
  # println("Error =$Error")
@test abs(Error[1]) < 1.0e-3
end
end
using .mmmPoissonmmmt10m2m
mmmPoissonmmmt10m2m.test()


module mPoisson_FE_Q4toT3_2
using FinEtools
using FinEtools.MeshExportModule
using Compat.Test
import LinearAlgebra: cholfact
function test()
  # println("""
  #
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit square, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular four-node QUADRILATERALS,
  # in a grid of 1000 x 1000 edges (1M quads, 1M degrees of freedom).
  # Version: 05/29/2017
  # """
  # )
  t0 = time()

  A = 1.0
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = -6.0; #heat source
  end
  tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
  N = 10;

  # println("Mesh generation")
  fens,fes = Q4block(A, A, N, N)
    fens,fes = Q4toT3(fens,fes)

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))


  # println("Searching nodes  for BC")
  l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
  l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
  l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
  l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
  List = vcat(l1, l2, l3, l4);
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)

  numberdofs!(Temp)

  t1 = time()

  m = MatHeatDiff(thermal_conductivity)
  femm = FEMMHeatDiff(IntegData(fes, TriRule(13)), m)

  # println("Conductivity")
  K=conductivity(femm, geom, Temp)
  #Profile.print()

  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  fi = ForceIntensity(FFlt, 1, getsource!);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")
  K = cholfact(K)
  # println("Solution of the factorized system")
  U=  K\(F1+F2)
  scattersysvec!(Temp, U[:])


  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  # using .meshExportModule

  # File =  "a.vtk"
  # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")

  Error = 0.0
  for k=1:size(fens.xyz,1)
    Error = Error .+ abs.(Temp.values[k,1] .- tempf(reshape(fens.xyz[k,:], (1,2))))
  end
  # println("Error =$Error")
  @test Error[1]<1.e-5

    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh(File, fes.conn, hcat(geom.values, Temp.values),
    # MeshExportModule.T3;
    # scalars=("Temperature", Temp.values))

  true
end
end
using .mPoisson_FE_Q4toT3_2
mPoisson_FE_Q4toT3_2.test()


module mPoisson_FE_T3_orientations
using FinEtools
using FinEtools.MeshExportModule
using Compat.Test
import LinearAlgebra: cholfact
function test()
    for orientation in [:a :b]

        # println("""
        #
        # Heat conduction example described by Amuthan A. Ramabathiran
        # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
        # Unit square, with known temperature distribution along the boundary,
        # and uniform heat generation rate inside.  Mesh of regular four-node QUADRILATERALS,
        # in a grid of 1000 x 1000 edges (1M quads, 1M degrees of freedom).
        # Version: 05/29/2017
        # """
        # )
        t0 = time()

        A = 1.0
        thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
        function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            forceout[1] = -6.0; #heat source
        end
        tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
        N = 10;

        # println("Mesh generation")
        fens,fes = T3block(A, A, N, N, orientation)

        geom = NodalField(fens.xyz)
        Temp = NodalField(zeros(size(fens.xyz,1),1))


        # println("Searching nodes  for BC")
        l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
        l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
        l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
        l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
        List = vcat(l1, l2, l3, l4);
        setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
        applyebc!(Temp)

        numberdofs!(Temp)

        t1 = time()

        m = MatHeatDiff(thermal_conductivity)
        femm = FEMMHeatDiff(IntegData(fes, TriRule(13)), m)

        # println("Conductivity")
        K=conductivity(femm, geom, Temp)
        #Profile.print()

        # println("Nonzero EBC")
        F2 = nzebcloadsconductivity(femm, geom, Temp);
        # println("Internal heat generation")
        fi = ForceIntensity(FFlt, 1, getsource!);
        F1 = distribloads(femm, geom, Temp, fi, 3);

        # println("Factorization")
        K = cholfact(K)
        # println("Solution of the factorized system")
        U=  K\(F1+F2)
        scattersysvec!(Temp, U[:])


        # println("Total time elapsed = $(time() - t0) [s]")
        # println("Solution time elapsed = $(time() - t1) [s]")

        # using .meshExportModule

        # File =  "a.vtk"
        # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")

        Error = 0.0
        for k=1:size(fens.xyz,1)
            Error = Error .+ abs.(Temp.values[k,1] .- tempf(reshape(fens.xyz[k,:], (1,2))))
        end
        # println("Error =$Error")
        @test Error[1]<1.e-5

        # File =  "a.vtk"
        # MeshExportModule.vtkexportmesh(File, fes.conn, hcat(geom.values, Temp.values),
        # MeshExportModule.T3;
        # scalars=("Temperature", Temp.values))

    end
end
end
using .mPoisson_FE_T3_orientations
mPoisson_FE_T3_orientations.test()

module mmT129b_l2_uqm
using FinEtools
using Compat.Test
import LinearAlgebra: cholfact
function test()
    L = 6.0;
    kappa = reshape([4.0], 1, 1);
    Q  = 0.1;
    n = 7; # number of elements
    crosssection = 1.0

    fens,fes = L2block(L, n)


    # Define boundary conditions
    l1  = selectnode(fens; box=[0. 0.], inflate = L/n/100.0)
    l2  = selectnode(fens; box=[L L], inflate = L/n/100.0)

    essential1 = FDataDict("node_list"=>vcat(l1, l2), "temperature"=> 0.0);
    material = MatHeatDiff(kappa)
    femm = FEMMHeatDiff(IntegData(fes, GaussRule(1, 2), crosssection), material)

    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
    setebc!(Temp, vcat(l1, l2), true, 1, 0.0)
    applyebc!(Temp)
    numberdofs!(Temp)

    K = conductivity(femm, geom, Temp)

    F2 = nzebcloadsconductivity(femm, geom, Temp);

    fi = ForceIntensity(FFlt[Q]);
    F1 = distribloads(femm, geom, Temp, fi, 3);

    K = cholfact(K)
    U = K\(F1+F2)
    scattersysvec!(Temp,U[:])

    # println("maximum(U)-0.1102 = $(maximum(U)-0.1102)")
@test abs(maximum(U)-0.1102) < 1.0e-4
end
end
using .mmT129b_l2_uqm
mmT129b_l2_uqm.test()

module mmannulusconv1
using FinEtools
using FinEtools.AlgoHeatDiffModule
using Compat.Test
function test()


    # println("""
    # Annular region, ingoing and outgoing flux. Temperature at one node prescribed.
    # Minimum/maximum temperature ~(+/-)0.500.
    # Mesh of serendipity quadrilaterals.
    # This version uses the FinEtools algorithm module.
    # Version: 05/29/2017
    # """)

    t0 = time()

    kappa = 0.2*[1.0 0; 0 1.0]; # conductivity matrix
    te1 = -0.5
    te2 = 0.6
    hconv1, hconv2 = 1000.0, 1000.0
    rin =  1.0;#internal radius
    rex =  2.0; #external radius
    nr = 7; nc = 90;
    Angle = 2*pi;
    thickness =  1.0;
    tolerance = min(rin/nr,  rin/nc/2/pi)/10000;

    fens, fes = Q4annulus(rin, rex, nr, nc, Angle)
    fens, fes = mergenodes(fens,  fes,  tolerance);
    edge_fes = meshboundary(fes);

    # The convection boundary condition is applied at two pieces of surface
    # Side 1
    l1 = selectelem(fens, edge_fes, box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
    el1femm = FEMMHeatDiffSurf(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)), hconv1)
    cbc1 = FDataDict("femm"=>el1femm, "ambient_temperature"=>te1)
    # Side 2
    l2=selectelem(fens,edge_fes,box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
    el2femm = FEMMHeatDiffSurf(IntegData(subset(edge_fes, l2),  GaussRule(1, 2)), hconv2)
    cbc2 = FDataDict("femm"=>el2femm, "ambient_temperature"=>te2)

    material = MatHeatDiff(kappa)
    femm = FEMMHeatDiff(IntegData(fes,  GaussRule(2, 2)),  material)
    region1 = FDataDict("femm"=>femm)

    # Make model data
    modeldata = FDataDict("fens"=>fens,
                     "regions"=>[region1],
                     "convection_bcs"=>[cbc1, cbc2]);

    # Call the solver
    modeldata = FinEtools.AlgoHeatDiffModule.steadystate(modeldata)
    geom=modeldata["geom"]
    Temp=modeldata["temp"]
    # println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")
    #
    # println("Total time elapsed = ",time() - t0,"s")
@test abs(minimum(Temp.values) - -0.5000094018945275) < 1.e-4
@test abs(maximum(Temp.values) - 0.600010400530557) < 1.e-4
    # # Postprocessing
    # File = "annulusmod.vtk"
    # vtkexportmesh(File, fes.conn, [geom.values Temp.values],
    # FinEtools.MeshExportModule.Q4; scalars=[("Temperature", Temp.values)])
    # @async run(`"paraview.exe" $File`)

end
end
using .mmannulusconv1
mmannulusconv1.test()
