module mmmmmmmmmmmm
using FinEtools
using Compat.Test
function test()
    println("""
    pnpConcreteColumnsymb_conv.
    Two-dimensional heat transfer with convection.
    """
    )
    
    a=2.5; dy=a/2*sin(15/180*pi); dx=a/2*cos(15/180*pi); Q=4.5; k=1.8; Dz=1.0;
    h= 5.;
    
    m = MatHeatDiff([k 0; 0 k])
    
    modeldata = nothing
    
    fens = FENodeSet([0 0; dx -dy; dx dy; 2*dx -2*dy; 2*dx 2*dy])
    fes = FESetT3([1 2 3; 2 4 5; 2 5 3]);
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
    
    applyebc!(Temp)
    numberdofs!(Temp)
    Temp.dofnums[1] = 3
    Temp.dofnums[2] = 1
    Temp.dofnums[3] = 2
    Temp.dofnums[4] = 5
    Temp.dofnums[5] = 4
    
    bfes = FESetL2([4 5]);
    
    cfemm = FEMMHeatDiffSurf(IntegData(bfes, GaussRule(1, 2), Dz), h)
    femm = FEMMHeatDiff(IntegData(fes, TriRule(1), Dz), m)
    fi = ForceIntensity(FFlt[Q]);
    F1 = distribloads(femm, geom, Temp, fi, 3);
    K = conductivity(femm, geom, Temp)
    F2 = nzebcloadsconductivity(femm, geom, Temp);
    H = surfacetransfer(cfemm, geom, Temp);
    show(H)
    # F3 = surfacetransferloads(femm, geom, temp, amb);
    Factor = cholesky(K+H)
    U = Factor\(F1+F2)
    scattersysvec!(Temp, U[:])
    display(Temp)
    
end
end
using .mmmmmmmmmmmm
mmmmmmmmmmmm.test()




