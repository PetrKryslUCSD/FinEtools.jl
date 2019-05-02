

# module mmdivmat2
# using FinEtools
# using FinEtools.FEMMDeforLinearBaseModule: infsup_gh
# using Test
# import LinearAlgebra: norm, cholesky
# function test()
# 	Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = ( 1.0, 1.0, 1.0, 1, 1, 1, :a)
# 	Ea, nua, alphaa = ( 1.0, 0.3, 0.0)
# 	fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt,
# 	   nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
# 	MR  =  DeforModelRed3D

# 	# Property and material
# 	material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

# 	femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)

# 	geom = NodalField(fens.xyz)
# 	u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
# 	numberdofs!(u)

# 	infsup_gh(femm, geom, u);
# end
# end
# using .mmdivmat2
# mmdivmat2.test()
