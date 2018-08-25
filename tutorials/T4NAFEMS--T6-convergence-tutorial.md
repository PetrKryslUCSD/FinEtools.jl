# Two-dimensional heat transfer with convection: convergence study

## Description

Consider a plate of uniform thickness, measuring 0.6 m by 1.0 m. On one
short edge the temperature is fixed at 100 °C, and on one long edge the
plate is perfectly insulated so that the heat flux is zero through that
edge. The other two edges are losing heat via convection to an ambient
temperature of 0 °C. The thermal conductivity of the plate is 52.0 W/(m
.°K), and the convective heat transfer coefficient is 750 W/(m^2.°K).
There is no internal generation of heat. Calculate the temperature 0.2 m
along the un-insulated long side, measured from the intersection with the
fixed temperature side (we call it point A below). The reference result is 18.25 °C.

The reference temperature at the point A  is 18.25 °C according to the
NAFEMS publication (which cites the book Carslaw, H.S. and J.C. Jaeger,
Conduction of Heat in Solids. 1959: Oxford University Press).

The present  tutorial will investigate the reference temperature  and it
will attempt to  estimate the  limit value more precisely using a
sequence of meshes and Richardson's extrapolation.

```julia
using FinEtools
using PGFPlotsX

kappa = [52. 0; 0 52.]*phun("W/(M*K)"); # conductivity matrix
h = 750*phun("W/(M^2*K)");# surface heat transfer coefficient
Width = 0.6*phun("M");# Geometrical dimensions
Height = 1.0*phun("M");
HeightA = 0.2*phun("M");
Thickness = 0.1*phun("M");
tolerance  = Width/1000;

material = MatHeatDiff(kappa);
```

The extrapolation  requires solutions  on a sequence of meshes. The model will be constructed repeatedly  inside the loop below.

```julia
resultsTempA = FFlt[]; params = FFlt[];
```

The model is constructed as follows: In order to  get a node at precisely the correct location  of point A, the domain is split into two rectangles which are meshed separately and then merged.  The meshes of the rectangles are then  repeatedly refined by quadri-section: each triangle is split into four triangles  of half the size.   The three node triangles are then converted  to  quadratic (six-node) triangles. The boundary of this mesh is then extracted and used  to construct the FEMM for the convection boundary condition. The model data is then passed to the steady-state solver. Finally, the temperature at point A is extracted and stored, and the  relative element size is also recorded.

```julia
let
    for nref = 3:5
```

The mesh is created from two rectangles, with two triangles each to begin with

```julia
        fens,fes = T3blockx([0.0, Width], [0.0, HeightA])
        fens2,fes2 = T3blockx([0.0, Width], [HeightA, Height])
        fens,newfes1,fes2 = mergemeshes(fens, fes, fens2, fes2, tolerance)
        fes = cat(newfes1,fes2)
```

Refine the mesh desired number of times

```julia
        for ref = 1:nref
            fens,fes = T3refine(fens,fes);
        end
```

Convert  the mesh of three-node triangles  to six-node triangles

```julia
        fens, fes = T3toT6(fens,fes);
```

Extract the boundary of the mesh of the interior

```julia
        bfes = meshboundary(fes)
```

Define boundary conditions:
The prescribed temperature is applied along edge 1 (the bottom
edge in Figure 1).

```julia
        l1 = selectnode(fens; box=[0. Width 0. 0.], inflate=tolerance)
        essential1 = FDataDict("node_list"=>l1, "temperature"=> 100.);
```

The convection boundary condition is applied along the edges
2,3,4. The elements along the boundary are quadratic line
elements L3. The order-four Gauss quadrature is sufficiently accurate.

```julia
        l2 = selectelem(fens, bfes; box=[Width Width  0.0 Height], inflate =tolerance)
        l3 = selectelem(fens, bfes; box=[0.0 Width Height Height], inflate =tolerance)
        cfemm = FEMMHeatDiffSurf(IntegData(subset(bfes,vcat(l2,l3)), GaussRule(1, 4), Thickness), h)
        convection1 = FDataDict("femm"=>cfemm, "ambient_temperature"=>0.);
```

FEMM  for the interior of the domain

```julia
        femm = FEMMHeatDiff(IntegData(fes, TriRule(3), Thickness), material)
        region1 = FDataDict("femm"=>femm)
```

Make the model data

```julia
        modeldata = FDataDict("fens"=> fens, "regions"=>[region1],
                                "essential_bcs"=>[essential1], "convection_bcs"=>[convection1]);
```

Call the solver

```julia
        modeldata = FinEtools.AlgoHeatDiffModule.steadystate(modeldata)
```

Find the node at the monitored location of the point A  [coordinates (Width,HeightA)].

```julia
        l4 = selectnode(fens; box=[Width Width HeightA HeightA], inflate =tolerance)
```

Collect the temperature and store it

```julia
        Temp = modeldata["temp"]
        push!(resultsTempA, Temp.values[l4][1]); push!(params, 1.0/2^nref)
    end
end
```

These are the computed results for the temperature at point A:

```julia
println("$( resultsTempA  )")
```

Richardson extrapolation can now be performed on the series of  data points.

```julia
solnestim, beta, c, residual = FinEtools.AlgoBaseModule.richextrapol(resultsTempA, params)
println("Solution estimate = $(solnestim)")
println("Convergence rate estimate  = $(beta )")
```

Finally we plot the estimated true error and the approximate error  (as the differences of the successive solutions). The slope of the approximate error curve  is very close to the slope of the estimated true error, which is an indication that the extrapolation  had a chance to work.

```julia
@pgf p = LogLogAxis({
                xlabel = "Element size",
                ylabel = "Temperature Error",
                title = "",
                "legend pos={south east}",
                xmin = "1.0e-2",
                xmax = "1.0",
                ymin = 1.0e-3, ymax = 0.1
            },
            Plot({
                color = "red",
                mark  = "o"
            }, Table([:x => params, :y => abs.(resultsTempA .- solnestim)])), LegendEntry("estimated true error"),
            Plot({
                color = "red",
                mark  = "x"
            }, Table([:x => params[1:2], :y => abs.(diff(resultsTempA))])), LegendEntry("approximate error"),
            );
display(p)
PGFPlotsX.print_tex(p)#-
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

