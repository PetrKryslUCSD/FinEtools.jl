

# module mmcubevibrationESNICET4
# using FinEtools
# using Test
# using Arpack: eigs
# import LinearAlgebra: norm, cholesky, cross
# function test()
# 	# println("""
# 	#         % Vibration modes of unit cube  of almost incompressible material.
# 	#         %
# 	#         % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
# 	#         % tetrahedral. International Journal for Numerical Methods in
# 	#         % Engineering 67: 841-867.""")
# 	#         t0 = time()


# 	E = 1*phun("PA");
# 	nu = 0.499;
# 	rho = 1*phun("KG/M^3");
# 	a = 1*phun("M"); b = a; h =  a;
# 	n1 = 24;# How many element edges per side?
# 	na =  n1; nb =  n1; nh  = n1;
# 	neigvs = 60                   # how many eigenvalues
# 	OmegaShift = (0.01*2*pi)^2;

# 	MR = DeforModelRed3D
# 	fens,fes  = T4block(a,b,h, na,nb,nh)
# 	# fens,fes  = T4refine20(fens,fes)

# 	geom = NodalField(fens.xyz)
# 	u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

# 	numberdofs!(u)

# 	material=MatDeforElastIso(MR, rho, E, nu, 0.0)

# 	femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
# 	femm = associategeometry!(femm, geom)

# 	K =stiffness(femm, geom, u)
# 	M =mass(femm, geom, u)
# 	d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
# 	d = d .- OmegaShift;
# 	fs = real(sqrt.(complex(d)))/(2*pi)
# 	println("Eigenvalues: $fs [Hz]")

# 	for mode = 7:neigvs
# 		scattersysvec!(u, v[:, mode])
# 		fld = fieldfromintegpoints(femm, geom, u, :pressure, 1; nodevalmethod = :averaging)
# 		File =  "cube-mode$(mode).vtk"
# 		vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.T4; vectors = [("u", u.values)], scalars = [("pressure", fld.values)])
# 	end
	
# 	@async run(`"paraview.exe" $File`)

# 	# @test abs(fs[7]-0.26259869196259) < 1.0e-5
# end

# end
# using .mmcubevibrationESNICET4
# mmcubevibrationESNICET4.test()


# module mmcubevibrationESNICEH8
# using FinEtools
# using Test
# using Arpack: eigs
# import LinearAlgebra: norm, cholesky, cross
# function test()
# 	# println("""
# 	#         % Vibration modes of unit cube  of almost incompressible material.
# 	#         %
# 	#         % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
# 	#         % tetrahedral. International Journal for Numerical Methods in
# 	#         % Engineering 67: 841-867.""")
# 	#         t0 = time()


# 	E = 1*phun("PA");
# 	nu = 0.499;
# 	rho = 1*phun("KG/M^3");
# 	a = 1*phun("M"); b = a; h =  a;
# 	n1 = 24;# How many element edges per side?
# 	na =  n1; nb =  n1; nh  = n1;
# 	neigvs = 20                   # how many eigenvalues
# 	OmegaShift = (0.01*2*pi)^2;

# 	MR = DeforModelRed3D
# 	fens,fes  = H8block(a,b,h, na,nb,nh)
# 	# fens,fes  = T4refine20(fens,fes)

# 	geom = NodalField(fens.xyz)
# 	u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

# 	numberdofs!(u)

# 	material=MatDeforElastIso(MR, rho, E, nu, 0.0)

# 	femm = FEMMDeforLinearESNICEH8(MR, IntegDomain(fes, NodalTensorProductRule(3)), material)
# 	femm = associategeometry!(femm, geom)

# 	K =stiffness(femm, geom, u)
# 	M =mass(femm, geom, u)
# 	d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
# 	d = d .- OmegaShift;
# 	fs = real(sqrt.(complex(d)))/(2*pi)
# 	println("Eigenvalues: $fs [Hz]")

# 	for mode = 7:neigvs
# 		scattersysvec!(u, v[:, mode])
# 		fld = fieldfromintegpoints(femm, geom, u, :pressure, 1; nodevalmethod = :averaging)
# 		File =  "cube-mode$(mode).vtk"
# 		vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.H8; vectors = [("u", u.values)], scalars = [("pressure", fld.values)])
# 	end
	
# 	@async run(`"paraview.exe" $File`)

# 	# @test abs(fs[7]-0.26259869196259) < 1.0e-5
# end

# end
# using .mmcubevibrationESNICEH8
# mmcubevibrationESNICEH8.test()


# module mmcubevibrationESNICEH8b
# using FinEtools
# using Test
# using Arpack: eigs
# import LinearAlgebra: norm, cholesky, cross
# function test()
# 	# println("""
# 	#         % Vibration modes of unit cube  of almost incompressible material.
# 	#         %
# 	#         % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
# 	#         % tetrahedral. International Journal for Numerical Methods in
# 	#         % Engineering 67: 841-867.""")
# 	#         t0 = time()


# 	E = 1*phun("PA");
# 	nu = 0.499;
# 	rho = 1*phun("KG/M^3");
# 	a = 1*phun("M"); b = a; h =  a;
# 	n1 = 9;# How many element edges per side?
# 	na =  n1; nb =  n1; nh  = n1;
# 	neigvs = 20                   # how many eigenvalues
# 	OmegaShift = (0.01*2*pi)^2;

# 	MR = DeforModelRed3D
# 	fens,fes  = T4block(a,b,h, na,nb,nh)
# 	fens,fes  = T4toH8(fens,fes)

# 	geom = NodalField(fens.xyz)
# 	u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

# 	numberdofs!(u)

# 	material=MatDeforElastIso(MR, rho, E, nu, 0.0)

# 	femm = FEMMDeforLinearESNICEH8(MR, IntegDomain(fes, NodalTensorProductRule(3)), material)
# 	femm = associategeometry!(femm, geom)

# 	K =stiffness(femm, geom, u)
# 	M =mass(femm, geom, u)
# 	d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
# 	d = d .- OmegaShift;
# 	fs = real(sqrt.(complex(d)))/(2*pi)
# 	println("Eigenvalues: $fs [Hz]")

# 	for mode = 7:neigvs
# 		scattersysvec!(u, v[:, mode])
# 		fld = elemfieldfromintegpoints(femm, geom, u, :pressure, 1; nodevalmethod = :averaging)
# 		File =  "cube-mode$(mode).vtk"
# 		vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.H8; vectors = [("u", u.values)], scalars = [("pressure", fld.values)])
# 	end
	
# 	@async run(`"paraview.exe" $File`)

# 	@test norm(fs.-[0.0, 0.0, 0.0, 1.76692e-8, 1.08952e-7, 1.40754e-7, 0.267052, 0.273879, 0.351477, 0.358597, 0.360482, 0.361265, 0.363118, 0.36379, 0.410238, 0.410963, 0.424434, 0.454963, 0.45869, 0.459053])./norm(fs) < 1.0e-5
# end

# end
# using .mmcubevibrationESNICEH8b
# mmcubevibrationESNICEH8b.test()

module minfsuptestik13
using FinEtools
using FinEtools.FEMMDeforLinearMSModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
using UnicodePlots
function test()
	lambdatol = sqrt(1e8*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [3, 4, 6]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
		fens, fes = T4toT10(fens, fes)
		
		for i = 1:count(fens)
			fens.xyz[i,:] = fens.xyz[i,:] + vec(reshape(fens.xyz[i,:], 1, 3)*A);
		end
		# @show fens.xyz
	
		# File =  "minfsuptest1.vtk"
		# vtkexportmesh(File, fens, fes)
		# try rm(File); catch end

		MR  =  DeforModelRed3D

		material = MatDeforElastIso(MR, E, nu)

		
		geom = NodalField(fens.xyz)
		u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
		bfes = meshboundary(fes)
		l1 = connectednodes(bfes)
		setebc!(u, l1, true, 1, 0.0)
		setebc!(u, l1, true, 2, 0.0)
		setebc!(u, l1, true, 3, 0.0)
		numberdofs!(u)

		femm  =  FEMMDeforLinearMST10(MR, IntegDomain(fes, TetRule(4)), material)
		femm = associategeometry!(femm, geom)
		Gh = infsup_gh(femm, geom, u);
		Sh = infsup_sh(femm, geom, u);

		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));

		# @show lambda
		abslambda = real.(filter(y -> !isnan(y), lambda));
		ix = findall(y  -> y < 0.0, abslambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));
		ix = findall(y  -> y > 0.0, abslambda);
		# pl = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		# display(pl)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	# @show lambdamin
	# pl = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	# display(pl)
	
	@test norm(lambdamin -  [0.0952635, 0.104529, 0.109738]) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptestik13
minfsuptestik13.test()