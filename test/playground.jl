
module minfsuptest1
using FinEtools
using FinEtools.FEMMDeforLinearBaseModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
# using UnicodePlots
function test()
	lambdatol = sqrt(1e4*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [2, 3, 4, 5, 6]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
		# fens, fes = T4toT10(fens, fes)
		# if parshiftmult != 0
		# 	ax=parshiftmult*Length;
		# 	ay=parshiftmult*Length;
		# 	az=parshiftmult*Length;
		# 	fens = transform_apply(fens,@(x, data) (x +[-ax/2*x(3)^2*x(2),ay/2*x(3)*x(1)^2,az/4*x(1)^3*x(2)+ax/4*x(3)^2*sin(0.2*x(2)^3)]), []);
		# end

		for i = 1:count(fens)
			fens.xyz[i,:] += A*fens.xyz[i,:];
		end

		# File =  "minfsuptest1.vtk"
		# vtkexportmesh(File, fens, fes)
		# try rm(File); catch end

		MR  =  DeforModelRed3D

		material = MatDeforElastIso(MR, E, nu)

		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(4)), material)

		geom = NodalField(fens.xyz)
		u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
		bfes = meshboundary(fes)
		l1 = connectednodes(bfes)
		setebc!(u, l1, true, 1, 0.0)
		setebc!(u, l1, true, 2, 0.0)
		setebc!(u, l1, true, 3, 0.0)
		numberdofs!(u)

		Gh = infsup_gh(femm, geom, u);
		Sh = infsup_sh(femm, geom, u);

		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));

		# @show lambda
		abslambda = deepcopy(lambda);
		ix = findall(y  -> y < 0.0, lambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, Length / ne)
	end

# @show lambdamin
# 	a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
# 	display(a)
	
	@test norm(lambdamin - [0.262065, 0.1709, 0.126159, 0.100228, 0.0828139]) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptest1
minfsuptest1.test()

