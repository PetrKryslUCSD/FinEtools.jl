# module mesh_rcm_1
# using FinEtools
# using FinEtools.MeshExportModule
# using FinEtools.MeshModificationModule: adjgraph, nodedegrees, revcm
# using Test
# using SparseArrays, UnicodePlots

# function test()
# 	conn = [9 1 8 4;
# 	1 3 2 8;
# 	8 2 7 5;
# 	2 6 7 7];
# 	nfens = 9;
# 	ag = adjgraph(conn, nfens)   
# 	nd = nodedegrees(ag)
# 	@test ag == Array{Int64,1}[[9, 8, 4, 3, 2], [1, 3, 8, 
# 	7, 5, 6], [1, 2, 8], [9, 1, 8], [8, 2, 7], [2, 7], [8, 2, 5, 6], [9, 1, 4, 3, 2, 7, 5], [1, 8, 4]]
# 	@test nd == [5, 6, 3, 3, 3, 2, 4, 7, 3] 
# 	numbering = revcm(ag, nd)
# 	@test numbering == [4, 9, 5, 8, 3, 1, 7, 2, 6]    

	
# 	A = sprand(nfens, nfens, 1/nfens)
# 	A = A+A'
# 	# display(spy(A))
# 	# display(spy(A[numbering, numbering]))
# end

# end
# using .mesh_rcm_1
# mesh_rcm_1.test()


# module mesh_rcm_2
# using FinEtools
# using FinEtools.MeshExportModule
# using FinEtools.MeshModificationModule: adjgraph, nodedegrees, revcm
# using Test
# using SparseArrays, UnicodePlots

# function test()
# 	nfens = 29;
	
# 	# A = sprand(nfens, nfens, 1/nfens)
# 	# A = A+A'
# 	A = spzeros(nfens, nfens)
# 	A[5 ,  1]  =  0.559559                                                
# 	A[24,  1]  =  0.079212                                                
# 	A[19,  2]  =  0.459102                                                
# 	A[8 ,  3]  =  0.844709                                                
# 	A[16,  3]  =  0.206808                                                
# 	A[24,  4]  =  0.82036                                                 
# 	A[1 ,  5]  =  0.559559                                                
# 	A[7 ,  5]  =  0.595562                                                
# 	A[28,  6]  =  0.151036                                                
# 	A[5 ,  7]  =  0.595562                                                
# 	A[3 ,  8]  =  0.844709                                                
# 	A[24,  8]  =  0.47093                                                 
# 	A[21,  9]  =  0.673707                                                
# 	A[22,  9]  =  0.492159                                                
# 	A[24,  9]  =  0.10736                                                 
# 	A[19, 10]  =  0.99992                                                 
# 	A[18, 11]  =  0.573561                                                
# 	A[22, 11]  =  0.803174                                                
# 	A[17, 12]  =  0.127183                                                
# 	A[28, 12]  =  0.644722                                                
# 	A[29, 14]  =  0.839357                                                
# 	A[3 , 16]  =  0.206808                                                
# 	A[20, 16]  =  0.0470789                                               
# 	A[12, 17]  =  0.127183                                                
# 	A[11, 18]  =  0.573561                                                
# 	A[2 , 19]  =  0.459102                                                
# 	A[10, 19]  =  0.99992                                                 
# 	A[16, 20]  =  0.0470789                                               
# 	A[26, 20]  =  0.661536                                                
# 	A[9 , 21]  =  0.673707                                                
# 	A[9 , 22]  =  0.492159                                                
# 	A[11, 22]  =  0.803174                                                
# 	A[24, 23]  =  0.373656                                                
# 	A[1 , 24]  =  0.079212                                                
# 	A[4 , 24]  =  0.82036                                                 
# 	A[8 , 24]  =  0.47093                                                 
# 	A[9 , 24]  =  0.10736                                                 
# 	A[23, 24]  =  0.373656                                                
# 	A[20, 26]  =  0.661536                                                
# 	A[6 , 28]  =  0.151036                                                
# 	A[12, 28]  =  0.644722                                                
# 	A[14, 29]  =  0.839357                                                

# 	ag = adjgraph(A)
# 	nd = nodedegrees(ag)
# 	numbering = revcm(ag, nd)
# 	display(spy(A))
# 	display(spy(A[numbering, numbering]))
# end

# end
# using .mesh_rcm_2
# mesh_rcm_2.test()


# module mesh_rcm_3
# using FinEtools
# using FinEtools.MeshExportModule
# using FinEtools.MeshModificationModule: adjgraph, nodedegrees, revcm
# using Test
# using SparseArrays, UnicodePlots

# function test()
# 	nfens = 19;
	
# 	# A1 = sprand(nfens, nfens, 1/nfens)
# 	# @show A1 = A1+A1'
# 	A1 = spzeros(nfens, nfens)
# 	A1[7 ,  1]  =  0.783374                                                
# 	A1[18,  1]  =  0.6411                                                  
# 	A1[8 ,  2]  =  0.66032                                                 
# 	A1[13,  2]  =  0.552169                                                
# 	A1[5 ,  3]  =  0.522678                                                
# 	A1[11,  4]  =  0.244274                                                
# 	A1[3 ,  5]  =  0.522678                                                
# 	A1[19,  5]  =  0.870687                                                
# 	A1[15,  6]  =  0.254443                                                
# 	A1[17,  6]  =  0.138423                                                
# 	A1[1 ,  7]  =  0.783374                                                
# 	A1[8 ,  7]  =  0.274651                                                
# 	A1[11,  7]  =  0.255421                                                
# 	A1[15,  7]  =  0.961861                                                
# 	A1[2 ,  8]  =  0.66032                                                 
# 	A1[7 ,  8]  =  0.274651                                                
# 	A1[11,  8]  =  0.0421145                                               
# 	A1[4 , 11]  =  0.244274                                                
# 	A1[7 , 11]  =  0.255421                                                
# 	A1[8 , 11]  =  0.0421145                                               
# 	A1[12, 11]  =  0.610131                                                
# 	A1[16, 11]  =  0.678996                                                
# 	A1[11, 12]  =  0.610131                                                
# 	A1[19, 12]  =  0.510702                                                
# 	A1[2 , 13]  =  0.552169                                                
# 	A1[18, 13]  =  0.0696182                                               
# 	A1[14, 14]  =  0.213021                                                
# 	A1[15, 14]  =  0.516788                                                
# 	A1[6 , 15]  =  0.254443                                                
# 	A1[7 , 15]  =  0.961861                                                
# 	A1[14, 15]  =  0.516788                                                
# 	A1[17, 15]  =  0.34131                                                 
# 	A1[11, 16]  =  0.678996                                                
# 	A1[6 , 17]  =  0.138423                                                
# 	A1[15, 17]  =  0.34131                                                 
# 	A1[1 , 18]  =  0.6411                                                  
# 	A1[13, 18]  =  0.0696182                                               
# 	A1[5 , 19]  =  0.870687                                                
# 	A1[12, 19]  =  0.510702
	
# 	A = vcat(hcat(A1, 0*A1, 0*A1), hcat(0*A1, 0*A1, 0*A1), hcat(0*A1, 0*A1, A1))
# 	ag = adjgraph(A)
# 	nd = nodedegrees(ag)
# 	numbering = revcm(ag, nd)
# 	display(spy(A))
# 	display(spy(A[numbering, numbering]))
# 	@test numbering == [55, 52, 44, 56, 51, 53, 39, 40, 45, 46, 54, 42, 49, 50, 57, 43, 41, 17, 14, 6, 18, 13, 15, 1, 2, 7, 8, 16, 4, 11, 12, 19, 5, 3, 48, 47, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 
# 26, 25, 24, 23, 22, 21, 20, 10, 9]
# end

# end
# using .mesh_rcm_3
# mesh_rcm_3.test()


# module mesh_rcm_5
# using FinEtools
# using FinEtools.MeshModificationModule: adjgraph, nodedegrees, revcm
# using Test
# using SparseArrays, UnicodePlots
# using LinearAlgebra
# using LinearAlgebra: norm

# function test()
# 	nfens = 511117;
	
# 	A = sprand(nfens, nfens, 1/nfens)
# 	A = A+A' + 1I
# 	@time ag = adjgraph(A)
# 	@time nd = nodedegrees(ag)
# 	@time numbering = revcm(ag, nd)
	
# 	# display(spy(A))
# 	# display(spy(A[numbering, numbering]))
# 	# b = rand(nfens)
# 	# inumbering = deepcopy(numbering)
# 	# inumbering[numbering] = 1:nfens
# 	# x = A * b
# 	# xp = A[numbering, numbering] * b[numbering]
# 	# xp[numbering]
# 	# @test norm(x - xp[inumbering]) < 1.0e-5*nfens
# end

# end
# using .mesh_rcm_5
# mesh_rcm_5.test()



