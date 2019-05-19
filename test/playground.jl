module minfsuptestik13
using FinEtools
using FinEtools.FEMMDeforLinearMSModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
# using UnicodePlots
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

		@show abslambda = sqrt.(sort(abslambda));
		# ix = findall(y  -> y > 0.0, abslambda);
		# pl = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		# display(pl)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	@show lambdamin
	# pl = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	# display(pl)
	
	@test norm(lambdamin -  [0.0952635, 0.104529, 0.109738]) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptestik13
minfsuptestik13.test()