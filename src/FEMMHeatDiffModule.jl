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
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.CSysModule: CSys, updatecsmat!
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!, gathervalues_asmat!
import FinEtools.NodalFieldModule: NodalField 
import FinEtools.ElementalFieldModule: ElementalField 
import FinEtools.AssemblyModule: SysvecAssemblerBase, SysmatAssemblerBase, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
import FinEtools.ForceIntensityModule: ForceIntensity
import FinEtools.FEMMBaseModule: FEMMAbstractBase, inspectintegpoints
import FinEtools.MatrixUtilityModule: add_gkgt_ut_only!, complete_lt!, locjac!
import LinearAlgebra: mul!, Transpose
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
import LinearAlgebra: norm, dot

# Type for heat diffusion finite element modeling machine.
mutable struct FEMMHeatDiff{S<:FESet, F<:Function, M<:MatHeatDiff} <: FEMMAbstractBase
    integdomain::IntegDomain{S, F} # geometry data
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
end

"""
    FEMMHeatDiff(integdomain::IntegDomain{S, F}, material::M) where {S<:FESet, F<:Function, M<:MatHeatDiff}

Construct with the default orientation matrix (identity).  
"""
function FEMMHeatDiff(integdomain::IntegDomain{S, F}, material::M) where {S<:FESet, F<:Function, M<:MatHeatDiff}
    return FEMMHeatDiff(integdomain, CSys(manifdim(integdomain.fes)), material)
end

function  buffers1(self::FEMMHeatDiff, geom::NodalField{FFlt}, temp::NodalField{FFlt})
    # Constants
    fes = self.integdomain.fes
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
    kappa_bar = fill(zero(FFlt), mdim, mdim); # buffer
    kappa_bargradNT = fill(zero(FFlt), mdim, nne); # buffer
    return dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat, elvec, elvecfix
end

"""
    conductivity(self::FEMMHeatDiff,
      assembler::A, geom::NodalField{FFlt},
      temp::NodalField{FFlt}) where {A<:SysmatAssemblerBase}

Compute the conductivity matrix.
"""
function conductivity(self::FEMMHeatDiff,  assembler::A, geom::NodalField{FFlt},  temp::NodalField{FFlt}) where {A<:SysmatAssemblerBase}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    # Prepare assembler and buffers
    dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat = buffers1(self, geom, temp)
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = self.material.tangentmoduli!(self.material, kappa_bar)
    startassembly!(assembler, size(elmat,1), size(elmat,2), count(fes), temp.nfreedofs, temp.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        fill!(elmat,  0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
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
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    # Prepare assembler and buffers
    dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat, elvec, elvecfix = buffers1(self, geom, temp)
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = self.material.tangentmoduli!(self.material, kappa_bar)
    startassembly!(assembler,  temp.nfreedofs);
    # Now loop over all finite elements in the set
    for i = 1:count(fes) # Loop over elements
        gatherfixedvalues_asvec!(temp, elvecfix, fes.conn[i]);# retrieve element coordinates
        if norm(elvecfix, Inf) != 0. # Is the load nonzero?
            fill!(elmat,  0.0);
            for j=1:npts # Loop over quadrature points
                locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
                Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
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

"""
    energy(self::FEMMHeatDiff, geom::NodalField{FFlt},  temp::NodalField{FFlt}) 

Compute the "energy" integral over the interior domain.
"""
function energy(self::FEMMHeatDiff, geom::NodalField{FFlt},  temp::NodalField{FFlt}) 
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    # Prepare assembler and buffers
    dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat, elvec, elvecfix = buffers1(self, geom, temp)
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar = self.material.tangentmoduli!(self.material, kappa_bar)
    gradT = fill(0.0, 1, size(gradN, 2))
    fluxT = deepcopy(gradT)
    energy = 0.0
    # Now loop over all finite elements in the set
    for i = 1:count(fes) # Loop over elements
        gathervalues_asvec!(temp, elvec, fes.conn[i]);# retrieve element coordinates
        for j=1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]);
            At_mul_B!(RmTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], RmTJ);
            At_mul_B!(gradT, elvec, gradN) 
            A_mul_B!(fluxT, gradT, kappa_bar) 
            energy += dot(vec(gradT), vec(fluxT)) * (Jac*w[j])
        end # Loop over quadrature points
    end
    return energy;
end

"""
    inspectintegpoints(self::FEMMHeatDiff, geom::NodalField{FFlt}, u::NodalField{T}, temp::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:heatflux; context...) where {F<:Function}

Inspect integration point quantities.

# Inputs
- `geom` - reference geometry field
- `u` - displacement field (ignored)
- `temp` - temperature field
- `felist` - indexes of the finite elements that are to be inspected: The fes
  to be included are: `fes[felist]`.
- `context`    - structure: see the update!() method of the material.
- `inspector` - function with the signature
        `idat = inspector(idat, j, conn, x, out, loc);`
   where
    `idat` - a structure or an array that the inspector may
           use to maintain some state,  for instance minimum or maximum of
           stress, `j` is the element number, `conn` is the element connectivity,
           `out` is the output of the update!() method,  `loc` is the location
           of the integration point in the *reference* configuration.
# Return
The updated inspector data is returned.
"""
function inspectintegpoints(self::FEMMHeatDiff, geom::NodalField{FFlt}, u::NodalField{T}, temp::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:heatflux; context...) where {T<:Number, F<:Function}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, RmTJ, gradN, kappa_bar, kappa_bargradNT, elmat, elvec, elvecfix = buffers1(self, geom, temp)
    # Thermal conductivity matrix is in local  material coordinates.
    kappa_bar =  self.material.thermal_conductivity;
    # Sort out  the output requirements
    outputcsys = self.mcsys; # default: report the vector quantities in the material coord system
    for apair in pairs(context)
        sy, val = apair
        if sy == :outputcsys
            outputcsys = val
        end
    end
    t= 0.0
    dt = 0.0
    Te = fill(zero(FFlt), nodesperelem(fes)) # nodal temperatures -- buffer
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    xe = fill(zero(FFlt), nne, sdim); # array of node coordinates -- buffer
    qpgradT = fill(zero(FFlt), 1, sdim); # Temperature gradient -- buffer
    qpflux = fill(zero(FFlt), sdim); # thermal strain -- buffer
    out1 = fill(zero(FFlt), sdim); # output -- buffer
    out =  fill(zero(FFlt), sdim);# output -- buffer
    # Loop over  all the elements and all the quadrature points within them
    for ilist = 1:length(felist) # Loop over elements
        i = felist[ilist];
        gathervalues_asmat!(geom, xe, fes.conn[i]);# retrieve element coords
        gathervalues_asvec!(temp, Te, fes.conn[i]);# retrieve element temperatures
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]);
            At_mul_B!(RmTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], RmTJ);
            # Quadrature point quantities
            A_mul_B!(qpgradT, reshape(Te, 1, :), gradN); # temperature gradient in material coordinates
            out = self.material.update!(self.material, qpflux, out, vec(qpgradT), 0.0, 0.0, loc, fes.label[i], quantity)
            if (quantity == :heatflux)   # Transform heat flux vector,  if that is "out"
                A_mul_B!(out1, transpose(self.mcsys.csmat), out);# To global coord sys
                A_mul_B!(out, outputcsys.csmat, out1);# To output coord sys
            end
            # Call the inspector
            idat = inspector(idat, i, fes.conn[i], xe, out, loc);
        end # Loop over quadrature points
    end # Loop over elements
    return idat; # return the updated inspector data
end

end
