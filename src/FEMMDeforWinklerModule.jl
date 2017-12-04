"""
    FEMMDeforWinklerModule

Module for operations on boundaries of domains to construct system matrices and
system vectors for linear deformation models with distributed-spring supports
(Winkler foundation model).
"""
module FEMMDeforWinklerModule

using FinEtools.FTypesModule
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: FESet, nodesperelem, manifdim
import FinEtools.IntegDataModule: IntegData, integrationdata, Jacobiansurface
import FinEtools.FieldModule: ndofs
import FinEtools.NodalFieldModule: NodalField 
import FinEtools.FEMMBaseModule: FEMMAbstractBase
import FinEtools.AssemblyModule: SysvecAssemblerBase, SysmatAssemblerBase, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!
import FinEtools.MatrixUtilityModule: add_nnt_ut_only!, complete_lt!, locjac!

"""
    FEMMDeforWinkler{S<:FESet, F<:Function}

Type for normal spring support  (Winkler).
"""
# Class for heat diffusion finite element modeling machine.
mutable struct FEMMDeforWinkler{S<:FESet, F<:Function} <: FEMMAbstractBase
    integdata::IntegData{S, F} # geometry data finite element modeling machine
end

"""
    surfacenormalspringstiffness(self::FEMMDeforWinkler, assembler::A,
                  geom::NodalField{FFlt}, u::NodalField{T},
                  springconstant::FFlt) where {A<:SysmatAssemblerBase, T<:Number}

Compute the stiffness matrix of surface normal spring.
Rationale: consider continuously distributed springs between the surface of the
solid body and the 'ground', in the direction normal to the surface. If the
spring coefficient becomes large, we have an approximate method of enforcing the
normal displacement to the surface.
"""
function surfacenormalspringstiffness(self::FEMMDeforWinkler, assembler::A,
    geom::NodalField{FFlt}, u::NodalField{T},
    springconstant::FFlt) where {A<:SysmatAssemblerBase, T<:Number}
    integdata = self.integdata
    # Constants
    nfes = count(integdata.fes); # number of finite elements in the set
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(integdata.fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(integdata.fes); # manifold dimension of the element
    Kedim = ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(integdata);
    # Prepare assembler and temporaries
    Ke = zeros(FFlt,Kedim,Kedim);                # element matrix -- used as a buffer
    dofnums = zeros(FInt, Kedim); # degree of freedom array -- used as a buffer
    loc = zeros(FFlt, 1,sdim); # quadrature point location -- used as a buffer
    J = fill(zero(FFlt), sdim,mdim); # Jacobian matrix -- used as a buffer
    startassembly!(assembler, Kedim, Kedim, nfes, u.nfreedofs, u.nfreedofs);
    for i = 1:nfes # Loop over elements
        fill!(Ke, 0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, integdata.fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobiansurface(integdata, J, loc, integdata.fes.conn[i], Ns[j]);
            n = surfacenormal(loc, J);# find the normal to the surface
            Nn = reshape(n*Ns[j]', Kedim, 1);# The normal n is a column vector
            add_nnt_ut_only!(Ke, Nn, springconstant*Jac*w[j])
        end # Loop over quadrature points
        complete_lt!(Ke)
        gatherdofnums!(u, dofnums, integdata.fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, Ke, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function surfacenormalspringstiffness(self::FEMMDeforWinkler,
              geom::NodalField{FFlt}, u::NodalField{T},
              springconstant::FFlt) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return surfacenormalspringstiffness(self, assembler, geom, u, springconstant);
end

"""
    surfacenormal(loc::FFltMat,J::FFltMat)

Compute local normal. This makes sense for bounding surfaces only.
"""
function  surfacenormal(loc::FFltMat, J::FFltMat)
    norml= zeros(FFlt, length(loc))
    # Produce a default normal
    if (size(J,1)==3) && (size(J,2)==2)# surface in three dimensions
        norml = cross(vec(J[:,1]),vec(J[:,2]));# outer normal to the surface
        norml = norml/norm(norml);
    elseif (size(J,1)==2)  && (size(J,2)==1)# curve in two dimensions
        norml= [J[2,1];-J[1,1]];# outer normal to the contour
        norml = norml/norm(norml);
    else
        error("No definition of normal vector");
    end
    return reshape(norml,length(norml),1) # return a column vector
end


end
