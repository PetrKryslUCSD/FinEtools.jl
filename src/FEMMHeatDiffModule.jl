"""
    FEMMHeatDiffModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear heat conduction/diffusion.
"""
module FEMMHeatDiffModule

export  FEMMHeatDiff
export conductivity,  nzebcloadsconductivity

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.FESetModule.gradN!
using FinEtools.CSysModule
using FinEtools.IntegDataModule
using FinEtools.FEMMBaseModule
using FinEtools.FieldModule
using FinEtools.NodalFieldModule
using FinEtools.ForceIntensityModule
using FinEtools.MatHeatDiffModule
using FinEtools.AssemblyModule
using FinEtools.MatrixUtilityModule.add_gkgt_ut_only!
using FinEtools.MatrixUtilityModule.complete_lt!
using FinEtools.MatrixUtilityModule.mv_product!
using FinEtools.MatrixUtilityModule: locjac!

# Type for heat diffusion finite element modeling machine.
mutable struct FEMMHeatDiff{S<:FESet, F<:Function, M<:MatHeatDiff} <: FEMMAbstractBase
    integdata::IntegData{S, F} # geometry data finite element modeling machine
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
end

"""
    FEMMHeatDiff(integdata::IntegData{S, F}, material::M) where {S<:FESet, F<:Function, M<:MatHeatDiff}

Construct with the default orientation matrix (identity).  
"""
function FEMMHeatDiff(integdata::IntegData{S, F}, material::M) where {S<:FESet, F<:Function, M<:MatHeatDiff}
    return FEMMHeatDiff(integdata, CSys(manifdim(integdata.fes)), material)
end

function  buffers(self::FEMMHeatDiff, geom::NodalField{FFlt}, temp::NodalField{FFlt})
    # Constants
    nfes = count(self.integdata.fes); # number of finite elements in the set
    ndn = ndofs(temp); # number of degrees of freedom per node
    nne = nodesperelem(self.integdata.fes); # number of nodes for element
    sdim = ndofs(geom);   # number of space dimensions
    mdim = manifdim(self.integdata.fes); # manifold dimension of the element
    Kedim = ndn*nne;      # dimension of the element matrix
    elmat = zeros(FFlt, Kedim, Kedim); # buffer
    elvec = zeros(FFlt, Kedim); # buffer
    elvecfix = zeros(FFlt, Kedim); # buffer
    conn = zeros(FInt, nne); # buffer
    x = zeros(FFlt, nne, sdim); # buffer
    dofnums = zeros(FInt, 1, Kedim); # buffer
    loc = zeros(FFlt, 1, sdim); # buffer
    J = eye(FFlt, sdim, mdim); # buffer
    RmTJ = zeros(FFlt, mdim, mdim); # buffer
    gradN = zeros(FFlt, nne, mdim); # buffer
    kappa_bargradNT = zeros(FFlt, mdim, nne); # buffer
    return conn, x, dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat, elvec, elvecfix
end

"""
    conductivity(self::FEMMHeatDiff,
      assembler::A, geom::NodalField{FFlt},
      temp::NodalField{FFlt}) where {A<:SysmatAssemblerBase}

Compute the conductivity matrix.
"""
function conductivity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:SysmatAssemblerBase}
    integdata = self.integdata
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(integdata);
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar =  self.material.thermal_conductivity;
    # Prepare assembler and buffers
    conn, x, dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat = buffers(self, geom, temp)
    startassembly!(assembler, size(elmat,1), size(elmat,2), count(integdata.fes), temp.nfreedofs, temp.nfreedofs);
    for i = 1:count(integdata.fes) # Loop over elements
        fill!(elmat,  0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, integdata.fes.conn[i], Ns[j], gradNparams[j])
            Jac = Jacobianvolume(integdata, J, loc, integdata.fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, integdata.fes.label[i]);
            At_mul_B!(RmTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
            gradN!(integdata.fes, gradN, gradNparams[j], RmTJ);
            # Add the product gradN*kappa_bar*gradNT*(Jac*w[j])
            add_gkgt_ut_only!(elmat, gradN, (Jac*w[j]), kappa_bar, kappa_bargradNT)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(temp, dofnums, integdata.fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function conductivity(self::FEMMHeatDiff, geom::NodalField{FFlt},  temp::NodalField{FFlt})
    assembler = SysmatAssemblerSparseSymm();
    return conductivity(self, assembler, geom, temp);
end

"""
    nzebcloadsconductivity(self::FEMMHeatDiff,
      assembler::A,  geom::NodalField{FFlt},
      temp::NodalField{FFlt}) where {A<:SysvecAssemblerBase}

Compute load vector for nonzero EBC of prescribed temperature.
"""
function nzebcloadsconductivity(self::FEMMHeatDiff, assembler::A,  geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:SysvecAssemblerBase}
    integdata = self.integdata
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(integdata);
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = self.material.thermal_conductivity;
    # Prepare assembler and buffers
    conn, x, dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat, elvec, elvecfix =
                buffers(self, geom, temp)
    startassembly!(assembler,  temp.nfreedofs);
    # Now loop over all finite elements in the set
    for i = 1:count(integdata.fes) # Loop over elements
        gatherfixedvalues_asvec!(temp, elvecfix, conn);# retrieve element coordinates
        if norm(elvecfix) != 0.     # Is the load nonzero?
            fill!(elmat,  0.0);
            for j=1:npts # Loop over quadrature points
                locjac!(loc, J, geom.values, integdata.fes.conn[i], Ns[j], gradNparams[j]) 
                Jac = Jacobianvolume(integdata, J, loc, integdata.fes.conn[i], Ns[j]);
                updatecsmat!(self.mcsys, loc, J, integdata.fes.label[i]);
                At_mul_B!(RmTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
                gradN!(integdata.fes, gradN, gradNparams[j], RmTJ);
                # Add the product gradN*kappa_bar*gradNT*(Jac*w[j])
                add_gkgt_ut_only!(elmat, gradN, (Jac*w[j]), kappa_bar, kappa_bargradNT)
            end # Loop over quadrature points
            complete_lt!(elmat)
            mv_product!(elvec, elmat, elvecfix) # compute  the load vector
            gatherdofnums!(temp, dofnums, integdata.fes.conn[i]); # retrieve degrees of freedom
            assemble!(assembler,  -elvec,  dofnums); # assemble element load vector
        end
    end
    return makevector!(assembler);
end

function nzebcloadsconductivity(self::FEMMHeatDiff,  geom::NodalField{FFlt},   temp::NodalField{FFlt})
    assembler = SysvecAssembler()
    return  nzebcloadsconductivity(self, assembler, geom, temp);
end


end
