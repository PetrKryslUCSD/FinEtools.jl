
module mPoisson_FE_T3_orientations
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
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
        thermal_conductivity = eye(2,2); # conductivity matrix
        function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            forceout[1] = -6.0; #heat source
        end
        tempf(x) = (1.0 + x[:,1].^2 + 2*x[:,2].^2);
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
        femm = FEMMHeatDiff(GeoD(fes, TriRule(13)), m)

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

        # using MeshExportModule

        # File =  "a.vtk"
        # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")

        Error = 0.0
        for k=1:size(fens.xyz,1)
            Error = Error+abs.(Temp.values[k,1]-tempf(reshape(fens.xyz[k,:], (1,2))))
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
using mPoisson_FE_T3_orientations
mPoisson_FE_T3_orientations.test()
