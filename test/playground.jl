
module minfsuptestik1
using FinEtools
using FinEtools.FEMMDeforLinearBaseModule: infsup_gh, infsup_sh
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
	for ne = [2, 3, 4]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
		fens, fes = T4toT10(fens, fes)
		# @show connasarray(fes)

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

		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)
		Gh = infsup_gh(femm, geom, u);
		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(4)), material)
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
		# a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		# display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	@show lambdamin
	# a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	# display(a)
	
	@test norm(lambdamin - [0.0729658, 0.0585958, 0.0459494] ) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptestik1
minfsuptestik1.test()


module minfsuptestik2
using FinEtools
using FinEtools.FEMMDeforLinearBaseModule: infsup_gh, infsup_sh
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
	for ne = [2, 3, 4]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
		# fens, fes = T4toT10(fens, fes)
		# @show connasarray(fes)

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

		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)
		Gh = infsup_gh(femm, geom, u);
		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)
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
		# a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		# display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	# @show lambdamin
	# a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	# display(a)
	
	@test norm(lambdamin - [0.270777, 0.179116, 0.132893]) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptestik2
minfsuptestik2.test()

module minfsuptestik3
using FinEtools
using FinEtools.FEMMDeforLinearBaseModule: infsup_gh, infsup_sh
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
	for ne = [2, 3, 4, 5, 6, 7, 8, 9]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = H8block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt)
		# @show connasarray(fes)

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

		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)
		Gh = infsup_gh(femm, geom, u);
		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)
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
		# a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		# display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	# @show lambdamin
	# a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	# display(a)
	
	@test norm(lambdamin - [0.447936, 0.317169, 0.305056, 0.300754, 0.291477, 0.290534, 
0.285866, 0.285624]) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptestik3
minfsuptestik3.test()

