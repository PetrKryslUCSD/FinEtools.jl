"""
    FEMMHeatDiffModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear heat conduction/diffusion.
"""
module FEMMHeatDiffModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: FESet, nodesperelem, manifdim, gradN!
import FinEtools.MatHeatDiffModule: MatHeatDiff
import FinEtools.IntegDataModule: IntegData, integrationdata, Jacobianvolume
import FinEtools.CSysModule: CSys, updatecsmat!
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!
import FinEtools.NodalFieldModule: NodalField 
import FinEtools.ElementalFieldModule: ElementalField 
import FinEtools.AssemblyModule: SysvecAssemblerBase, SysmatAssemblerBase, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
import FinEtools.ForceIntensityModule: ForceIntensity
import FinEtools.FEMMBaseModule: FEMMAbstractBase
import FinEtools.MatrixUtilityModule: add_gkgt_ut_only!, complete_lt!, locjac!

# Type for heat diffusion finite element modeling machine.
mutable struct FEMMHeatDiff{S<:FESet, F<:Function, M<:MatHeatDiff} <: FEMMAbstractBase
    integdata::IntegData{S, F} # geometry data
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
    fes = self.integdata.fes
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(temp); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);   # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    Kedim = ndn*nne;      # dimension of the element matrix
    elmat = fill(zero(FFlt), Kedim, Kedim); # buffer
    elvec = fill(zero(FFlt), Kedim); # buffer
    elvecfix = fill(zero(FFlt), Kedim); # buffer
    dofnums = fill(zero(FInt), Kedim); # buffer
    loc = fill(zero(FFlt), 1, sdim); # buffer
    J = fill(zero(FFlt), sdim, mdim); # buffer
    RmTJ = fill(zero(FFlt), mdim, mdim); # buffer
    gradN = fill(zero(FFlt), nne, mdim); # buffer
    kappa_bargradNT = fill(zero(FFlt), mdim, nne); # buffer
    return dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat, elvec, elvecfix
end

"""
    conductivity(self::FEMMHeatDiff,
      assembler::A, geom::NodalField{FFlt},
      temp::NodalField{FFlt}) where {A<:SysmatAssemblerBase}

Compute the conductivity matrix.
"""
function conductivity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:SysmatAssemblerBase}
    fes = self.integdata.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdata);
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar =  self.material.thermal_conductivity;
    # Prepare assembler and buffers
    dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat = buffers(self, geom, temp)
    startassembly!(assembler, size(elmat,1), size(elmat,2), count(fes), temp.nfreedofs, temp.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        fill!(elmat,  0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdata, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]);
            At_mul_B!(RmTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], RmTJ);
            # Add the product gradN*kappa_bar*gradNT*(Jac*w[j])
            add_gkgt_ut_only!(elmat, gradN, (Jac*w[j]), kappa_bar, kappa_bargradNT)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(temp, dofnums, fes.conn[i]);# retrieve degrees of freedom
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
    fes = self.integdata.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdata);
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = self.material.thermal_conductivity;
    # Prepare assembler and buffers
    dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat, elvec, elvecfix = buffers(self, geom, temp)
    startassembly!(assembler,  temp.nfreedofs);
    # Now loop over all finite elements in the set
    for i = 1:count(fes) # Loop over elements
        gatherfixedvalues_asvec!(temp, elvecfix, fes.conn[i]);# retrieve element coordinates
        if norm(elvecfix, Inf) != 0. # Is the load nonzero?
            fill!(elmat,  0.0);
            for j=1:npts # Loop over quadrature points
                locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
                Jac = Jacobianvolume(self.integdata, J, loc, fes.conn[i], Ns[j]);
                updatecsmat!(self.mcsys, loc, J, fes.label[i]);
                At_mul_B!(RmTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
                gradN!(fes, gradN, gradNparams[j], RmTJ);
                # Add the product gradN*kappa_bar*gradNT*(Jac*w[j])
                add_gkgt_ut_only!(elmat, gradN, (Jac*w[j]), kappa_bar, kappa_bargradNT)
            end # Loop over quadrature points
            complete_lt!(elmat)
            A_mul_B!(elvec, elmat, elvecfix) # compute  the load vector
            gatherdofnums!(temp, dofnums, fes.conn[i]); # retrieve degrees of freedom
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
