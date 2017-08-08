module mmmmmmmmmNAFEMSm
using FinEtools
using Base.Test
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
    cfemm = FEMMHeatDiffSurf(GeoD(subset(bfes,vcat(l2,l3)),
      GaussRule(1, 3), Thickness), h)
    convection1 = FDataDict("femm"=>cfemm, "ambient_temperature"=>0.);

    # The interior
    femm = FEMMHeatDiff(GeoD(fes, TriRule(3), Thickness), m)
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
  vtkexportmesh("T4NAFEMS--T6.vtk", regions[1]["femm"].geod.fes.conn,
  [geom.values Temp.values/100], FinEtools.MeshExportModule.T6;
  scalars=[("Temperature", Temp.values)])
  try  rm("T4NAFEMS--T6.vtk"); catch end
  vtkexportmesh("T4NAFEMS--T6--base.vtk", regions[1]["femm"].geod.fes.conn,
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
using mmmmmmmmmNAFEMSm
mmmmmmmmmNAFEMSm.test()
