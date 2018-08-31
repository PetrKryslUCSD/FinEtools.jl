module baffled_piston_examples
using FinEtools
using LinearAlgebra
using PGFPlotsX

function baffled_piston_H8_transient()
    rho = 1.21*phun("kg/m^3");# mass density
c  = 343.0*phun("m/s");# sound speed
bulk =  c^2*rho;
frequency=2000; # Hz
omega=2*pi*frequency;
P_piston = 1.0e-3*phun("MPa"); # amplitude of the piston pressure
R = 50.0*phun("mm");# radius of the piston
Ro = 400.0*phun("mm"); # radius of the enclosure
nref = 2;#number of refinements of the sphere around the piston
nlayers = 48;                     # number of layers of elements surrounding the piston
tolerance = R/(2^nref)/10
dt=1.0/frequency/20;
tfinal=50*dt;
nsteps=round(tfinal/dt)+1;

println("""

Baffled piston in a half-sphere domain with hard boundary condition.

Hexahedral mesh. Algorithm version.
""")

t0  =  time()

# Hexahedral mesh
fens,fes  =  H8sphere(R,nref);
bfes  =  meshboundary(fes)
# File  =   "baffledabc_boundary.vtk"
# vtkexportmesh(File, bfes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
#  @async run(`"paraview.exe" $File`)

l = selectelem(fens,bfes,facing = true,direction = [1.0 1.0  1.0])
ex(xyz, layer) = (R+layer/nlayers*(Ro-R))*xyz/norm(xyz)
fens1,fes1  =  H8extrudeQ4(fens, subset(bfes,l), nlayers, ex);
fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
fes = cat(newfes1,fes2)

# Piston surface mesh
bfes  =  meshboundary(fes)
l1 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])
l2 = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0], inflate = tolerance)
piston_fes = subset(bfes,intersect(l1,l2));
# File  =   "piston_fes.vtk"
# vtkexportmesh(File, piston_fes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
#  @async run(`"paraview.exe" $File`)

# Outer spherical boundary
louter = selectelem(fens, bfes, facing = true, direction = [1.0 1.0  1.0], dotmin= 0.001)
outer_fes = subset(bfes,louter);
# File  =   "outer_fes.vtk"
# vtkexportmesh(File, outer_fes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
#  @async run(`"paraview.exe" $File`)

println("Pre-processing time elapsed  =  ",time() - t0,"s")

t1  =  time()

material = MatAcoustFluid(bulk, rho)
# Region of the fluid
femm = FEMMAcoust(IntegData(fes, GaussRule(3, 2)), material)

geom = NodalField(fens.xyz);
P = NodalField(zeros(nnodes(geom),1));
piston_fenids = connectednodes(piston_fes);
setebc!(P, piston_fenids, true, 1, 0.0);
applyebc!(P);
numberdofs!(P);

S  =  acousticstiffness(femm, geom, P);
C  =  acousticmass(femm, geom, P);
D= (2.0/dt)*S +(dt/2.0)*C;

P0 = deepcopy(P)
Pdd0 = deepcopy(P)
P1 = deepcopy(P)
Pdd1 = deepcopy(P)
TMPF = deepcopy(P)
P = nothing # we don't need this field anymore
vP0 = zeros(P0.nfreedofs)
vP1 = deepcopy(vP0)
vPd0 = deepcopy(vP0)
vPd1 = deepcopy(vP0)

t = 0.0
P0.fixed_values[piston_fenids,1] .= P_piston*sin(omega*t)
Pdd0.fixed_values[piston_fenids,1] .= P_piston*(-omega^2)*sin(omega*t)
vP0 = gathersysvec!(P0, vP0)

nh = selectnode(fens, nearestto = [R+Ro/2, 0.0, 0.0] )
Pnh = [P1.values[nh, 1][1]]
ts = [t]
step =0;
while t <=tfinal
    step = step  +1;
    t=t+dt;
    P1.fixed_values[piston_fenids,1] .= P_piston*sin(omega*t)
    Pdd1.fixed_values[piston_fenids,1] .= P_piston*(-omega^2)*sin(omega*t)
    TMPF.fixed_values = P0.fixed_values + P1.fixed_values
    F = nzebcloadsacousticmass(femm, geom, TMPF);
    TMPF.fixed_values = Pdd0.fixed_values + Pdd1.fixed_values
    F = F + nzebcloadsacousticstiffness(femm, geom, TMPF);
    # println("$(norm(F))")
    vPd1 = D\((2/dt)*(S*vPd0) - C*(2*vP0+(dt/2)*vPd0) + F);
    vP1 = vP0 + (dt/2)*(vPd0+vPd1);
    scattersysvec!(P1, vP1); # store current pressure
    # Swap variables for the next step
    copyto!(vP0, vP1)
    copyto!(vPd0, vPd1)
    copyto!(P0, P1)
    copyto!(Pdd0, Pdd1)
    # Graphics output
    File  =   "baffled_piston-$(step).vtk"
    vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.H8;
        scalars = [(  "P", P1.values)])
    #  @async run(`"paraview.exe" $File`)
    push!(Pnh, P1.values[nh, 1][1])
    push!(ts, t)
    println("step = $( step )")
end

@pgf a = Axis({
    xlabel = "Time",
    ylabel = "Pnh",
    title = "Transient pressure"
},
Plot(Table([:x => vec(ts), :y => vec(Pnh)])))
display(a)

# using Plots
# plotly()
# plot(vec(ts), vec(Pnh))
# gui()
println("$(Pnh[end])")

# # Surface for the ABC
# abc1  =  FDataDict("femm"=>FEMMAcoustSurf(IntegData(outer_fes, GaussRule(2, 2)),
#           material))
#
# # Surface of the piston
# flux1  =  FDataDict("femm"=>FEMMAcoustSurf(IntegData(piston_fes, GaussRule(2, 2)),
#           material),  "normal_flux"=> -rho*a_piston+0.0im);
#
# # Make model data
# modeldata =  FDataDict("fens"=>  fens,
#                  "omega"=>omega,
#                  "regions"=>[region1],
#                  "flux_bcs"=>[flux1], "ABCs"=>[abc1])
#
# # Call the solver
# modeldata = FinEtools.AlgoAcoustModule.steadystate(modeldata)
#
# println("Computing time elapsed  =  ",time() - t1,"s")
# println("Total time elapsed  =  ",time() - t0,"s")
#
# geom = modeldata["geom"]
# P = modeldata["P"]
#
# File  =   "baffledabc.vtk"
# vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
# scalars = [("absP", abs.(P.values))])
# @async run(`"paraview.exe" $File`)
#
# # using Winston
# # pl  =  FramedPlot(title = "Matrix",xlabel = "x",ylabel = "Re P, Im P")
# # setattr(pl.frame, draw_grid = true)
# # add(pl, Curve([1:length(C[:])],vec(C[:]), color = "blue"))
#
# # # pl = plot(geom.values[nLx,1][ix],scalars[nLx][ix])
# # # xlabel("x")
# # # ylabel("Pressure")
# # display(pl)
#
# true
end # baffled_piston_H8_transient

function allrun()
    println("#####################################################") 
    println("# baffled_piston_H8_transient ")
    baffled_piston_H8_transient()
    return true
end # function allrun

end # module baffled_piston_examples
