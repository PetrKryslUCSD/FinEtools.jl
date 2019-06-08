"""
    AbstractFEMMDeforLinearBaseModule

Base module for operations on interiors of domains to construct system matrices and
system vectors for linear deformation models.
"""
module FEMMDeforLinearBaseModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: gradN!, nodesperelem, manifdim
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!, gathervalues_asmat!
import FinEtools.NodalFieldModule: NodalField, nnodes 
import FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
import FinEtools.FEMMBaseModule: AbstractFEMM, inspectintegpoints
import FinEtools.CSysModule: CSys, updatecsmat!
import FinEtools.DeforModelRedModule: nstressstrain, nthermstrain, Blmat!, divmat, vgradmat
import FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, add_btv!, locjac!, add_nnt_ut_only!
import FinEtools.MatDeforModule: rotstressvec!
import FinEtools.MatModule: massdensity
import FinEtools.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
import FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
import LinearAlgebra: Transpose, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
import LinearAlgebra: norm, dot

"""
    AbstractFEMMDeforLinear <: AbstractFEMMBase

Abstract type of FEMM for linear deformation.
"""
abstract type AbstractFEMMDeforLinear <: AbstractFEMM end

function _buffers(self::AbstractFEMMDeforLinear, geom::NodalField, u::NodalField)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare _buffers
    elmat = fill(zero(FFlt), elmatdim, elmatdim);      # element matrix -- buffer
    dofnums = zeros(FInt, elmatdim); # degree of freedom array -- buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    gradN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    D = fill(zero(FFlt), nstrs, nstrs); # material stiffness matrix -- buffer
    B = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    DB = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    elvecfix = fill(zero(FFlt), elmatdim); # vector of prescribed displ. -- buffer
    elvec = fill(zero(FFlt), elmatdim); # element vector -- buffer
    return dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix
end

"""
    mass(self::AbstractFEMMDeforLinear,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the consistent mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function mass(self::AbstractFEMMDeforLinear,  assembler::A,  geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat = _buffers(self, geom, u)  # Prepare buffers
    rho::FFlt = massdensity(self.material); # mass density
    NexpTNexp = FFltMat[];# basis f. matrix -- buffer
    ndn = ndofs(u)
    Indn = [i==j ? one(FFlt) : zero(FFlt) for i=1:ndn, j=1:ndn] # "identity"
    for j = 1:npts # This quantity is the same for all quadrature points
        Nexp = fill(zero(FFlt), ndn, size(elmat,1))
        for l1 = 1:nodesperelem(fes)
            Nexp[1:ndn, (l1-1)*ndn+1:(l1)*ndn] = Indn * Ns[j][l1];
        end
        push!(NexpTNexp, Nexp'*Nexp);
    end
    startassembly!(assembler,  size(elmat,1),  size(elmat,2),  count(fes), u.nfreedofs,  u.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        fill!(elmat,  0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            thefactor::FFlt =(rho*Jac*w[j]);
            elmat .+= NexpTNexp[j]*thefactor
        end # Loop over quadrature points
        gatherdofnums!(u,  dofnums,  fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler,  elmat,  dofnums,  dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function mass(self::AbstractFEMMDeforLinear,  geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return mass(self, assembler, geom, u);
end

"""
    stiffness(self::AbstractFEMMDeforLinear, assembler::A,
          geom::NodalField{FFlt},
          u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute and assemble  stiffness matrix.
"""
function stiffness(self::AbstractFEMMDeforLinear, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix = _buffers(self, geom, u)
    tangentmoduli!(self.material, D, 0.0, 0.0, loc, 0)
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), u.nfreedofs, u.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        fill!(elmat,  0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]);
            At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], csmatTJ);
            Blmat!(self.mr, B, Ns[j], gradN, loc, self.mcsys.csmat);
            add_btdb_ut_only!(elmat, B, Jac*w[j], D, DB)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(u, dofnums, fes.conn[i]); # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function stiffness(self::AbstractFEMMDeforLinear, geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom, u);
end

"""
    nzebcloadsstiffness(self::AbstractFEMMDeforLinear,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysvecAssembler, T<:Number}

Compute load vector for nonzero EBC for fixed displacement.
"""
function nzebcloadsstiffness(self::AbstractFEMMDeforLinear,  assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysvecAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix = _buffers(self, geom, u) 
    tangentmoduli!(self.material, D, 0.0, 0.0, loc, 0)
    startassembly!(assembler,  u.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gatherfixedvalues_asvec!(u, elvecfix, fes.conn[i]);# retrieve element displacement vector
        if norm(elvecfix, Inf) != 0.0   # Is the load nonzero?
            fill!(elmat,  0.0); # Initialize element matrix
            for j = 1:npts # Loop over quadrature points
                locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
                Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
                updatecsmat!(self.mcsys, loc, J, fes.label[i]);
                At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
                gradN!(fes, gradN, gradNparams[j], csmatTJ);
                Blmat!(self.mr, B, Ns[j], gradN, loc, self.mcsys.csmat);
                add_btdb_ut_only!(elmat, B, Jac*w[j], D, DB)
            end # Loop over quadrature points
            complete_lt!(elmat)
            gatherdofnums!(u, dofnums, fes.conn[i]); # retrieve degrees of freedom
            A_mul_B!(elvec, elmat, elvecfix)
            assemble!(assembler,  -elvec,  dofnums); # assemble element load vector
        end
    end # Loop over elements
    return makevector!(assembler);
end

function nzebcloadsstiffness(self::AbstractFEMMDeforLinear, geom::NodalField{FFlt}, u::NodalField{T}) where {T<:Number}
    assembler = SysvecAssembler()
    return  nzebcloadsstiffness(self, assembler, geom, u);
end

"""
    thermalstrainloads(self::AbstractFEMMDeforLinear, assembler::A,
        geom::NodalField{FFlt}, u::NodalField{T},
        dT::NodalField{FFlt}) where {A<:AbstractSysvecAssembler, T<:Number}

Compute the thermal-strain load vector.
"""
function  thermalstrainloads(self::AbstractFEMMDeforLinear, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}, dT::NodalField{FFlt}) where {A<:AbstractSysvecAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix = _buffers(self, geom, u)
    t= 0.0
    dt = 0.0
    DeltaT = fill(zero(FFlt), nodesperelem(fes))
    strain = fill(zero(FFlt), nstressstrain(self.mr)); # total strain -- buffer
    thstrain = fill(zero(FFlt), nthermstrain(self.mr)); # thermal strain -- buffer
    thstress = fill(zero(FFlt), nstressstrain(self.mr)); # thermal stress -- buffer
    startassembly!(assembler,  u.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asvec!(dT, DeltaT, fes.conn[i]);# retrieve element temperatures
        if norm(DeltaT, Inf) != 0     # Is the thermal increment nonzero?
            fill!(elvec,  0.0); # Initialize element matrix
            for j = 1:npts # Loop over quadrature points
                locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
                Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
                updatecsmat!(self.mcsys, loc, J, fes.label[i]);
                At_mul_B!(csmatTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
                gradN!(fes, gradN, gradNparams[j], csmatTJ);#Do: gradN = gradNparams[j]/csmatTJ;
                Blmat!(self.mr, B, Ns[j], gradN, loc, self.mcsys.csmat);# strains in mater cs, displ in global cs
                thermalstrain!(self.material, thstrain, dot(vec(Ns[j]), DeltaT))
                thstress = update!(self.material, thstress, thstress, strain, thstrain, t, dt, loc, fes.label[i], :nothing)
                add_btv!(elvec, B, thstress, (-1)*(Jac*w[j]))
            end
            gatherdofnums!(u, dofnums, fes.conn[i]); # retrieve degrees of freedom
            assemble!(assembler,  elvec,  dofnums); # assemble element load vector
        end
    end # Loop over elements
    return makevector!(assembler);
end

function thermalstrainloads(self::AbstractFEMMDeforLinear, geom::NodalField{FFlt}, u::NodalField{T}, dT::NodalField{FFlt}) where {T<:Number}
    assembler = SysvecAssembler()
    return  thermalstrainloads(self, assembler, geom, u, dT);
end

"""
    inspectintegpoints(self::AbstractFEMMDeforLinear,
      geom::NodalField{FFlt},  u::NodalField{T},
      dT::NodalField{FFlt},
      felist::FIntVec,
      inspector::F,  idat, quantity=:Cauchy;
      context...) where {T<:Number, F<:Function}

Inspect integration point quantities.

- `geom` - reference geometry field
- `u` - displacement field
- `dT` - temperature difference field
- `felist` - indexes of the finite elements that are to be inspected:
     The fes to be included are: `fes[felist]`.
- `context`    - structure: see the update!() method of the material.
- `inspector` - functionwith the signature
        idat = inspector(idat, j, conn, x, out, loc);
   where
    `idat` - a structure or an array that the inspector may
           use to maintain some state,  for instance minimum or maximum of
           stress, `j` is the element number, `conn` is the element connectivity,
           `out` is the output of the update!() method,  `loc` is the location
           of the integration point in the *reference* configuration.
### Return
The updated inspector data is returned.
"""
function inspectintegpoints(self::FEMM, geom::NodalField{FFlt},  u::NodalField{T}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {FEMM<:AbstractFEMMDeforLinear, T<:Number, F<:Function}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix = _buffers(self, geom, u)
    # Sort out  the output requirements
    outputcsys = self.mcsys; # default: report the stresses in the material coord system
    for apair in pairs(context)
        sy, val = apair
        if sy == :outputcsys
            outputcsys = val
        end
    end
    t= 0.0
    dt = 0.0
    dTe = fill(zero(FFlt), nodesperelem(fes)) # nodal temperatures -- buffer
    ue = fill(zero(FFlt), size(elmat, 1)); # array of node displacements -- buffer
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    xe = fill(zero(FFlt), nne, sdim); # array of node coordinates -- buffer
    qpdT = 0.0; # node temperature increment
    qpstrain = fill(zero(FFlt), nstressstrain(self.mr), 1); # total strain -- buffer
    qpthstrain = fill(zero(FFlt), nthermstrain(self.mr)); # thermal strain -- buffer
    qpstress = fill(zero(FFlt), nstressstrain(self.mr)); # stress -- buffer
    out1 = fill(zero(FFlt), nstressstrain(self.mr)); # stress -- buffer
    out =  fill(zero(FFlt), nstressstrain(self.mr));# output -- buffer
    # Loop over  all the elements and all the quadrature points within them
    for ilist = 1:length(felist) # Loop over elements
        i = felist[ilist];
        gathervalues_asmat!(geom, xe, fes.conn[i]);# retrieve element coords
        gathervalues_asvec!(u, ue, fes.conn[i]);# retrieve element displacements
        gathervalues_asvec!(dT, dTe, fes.conn[i]);# retrieve element temp. increments
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]);
            At_mul_B!(csmatTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], csmatTJ);
            Blmat!(self.mr, B, Ns[j], gradN, loc, self.mcsys.csmat);
            updatecsmat!(outputcsys, loc, J, fes.label[i]);
            # Quadrature point quantities
            A_mul_B!(qpstrain, B, ue); # strain in material coordinates
            qpdT = dot(vec(dTe), vec(Ns[j]));# Quadrature point temperature increment
            thermalstrain!(self.material, qpthstrain, qpdT)
            # Material updates the state and returns the output
            out = update!(self.material, qpstress, out, vec(qpstrain), qpthstrain, t, dt, loc, fes.label[i], quantity)
            if (quantity == :Cauchy)   # Transform stress tensor,  if that is "out"
                (length(out1) >= length(out)) || (out1 = zeros(length(out)))
                rotstressvec!(self.mr, out1, out, transpose(self.mcsys.csmat))# To global coord sys
                rotstressvec!(self.mr, out, out1, outputcsys.csmat)# To output coord sys
            end
            # Call the inspector
            idat = inspector(idat, i, fes.conn[i], xe, out, loc);
        end # Loop over quadrature points
    end # Loop over elements
    return idat; # return the updated inspector data
end

function inspectintegpoints(self::FEMM, geom::NodalField{FFlt},  u::NodalField{T}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {FEMM<:AbstractFEMMDeforLinear, T<:Number, F<:Function}
    dT = NodalField(fill(zero(FFlt), nnodes(geom), 1)) # zero difference in temperature
    return inspectintegpoints(self, geom, u, dT, felist, inspector, idat, quantity; context...);
end

function _buffers2(self::AbstractFEMMDeforLinear, geom::NodalField, u::NodalField)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare buffers
    elmat = fill(zero(FFlt), elmatdim, elmatdim);      # element matrix -- buffer
    dofnums = zeros(FInt, elmatdim); # degree of freedom array -- buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    gradN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
   	return dofnums, loc, J, csmatTJ, gradN, elmat
end

"""
    infsup_gh(self::AbstractFEMMDeforLinear, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the matrix to produce the norm of the divergence of the displacement.

This matrix is used in the numerical infsup test (Klaus-Jurgen Bathe, The
inf-sup condition and its evaluation for mixed finite element methods,
Computers and Structures 79 (2001) 243-252.)

!!! note 
This computation has not been optimized in any way. It can be expected to be
inefficient.
"""
function infsup_gh(self::AbstractFEMMDeforLinear, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
	fes = self.integdomain.fes
	npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
	dofnums, loc, J, csmatTJ, gradN, elmat = _buffers2(self, geom, u)
	startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), u.nfreedofs, u.nfreedofs);
	for i = 1:count(fes) # Loop over elements
	    fill!(elmat,  0.0); # Initialize element matrix
	    for j = 1:npts # Loop over quadrature points
	        locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
	        Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
	        updatecsmat!(self.mcsys, loc, J, fes.label[i]);
	        At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
	        gradN!(fes, gradN, gradNparams[j], csmatTJ);
	        divm = divmat(self.mr, Ns[j], gradN, loc);
	        elmat += (transpose(divm) * divm) * (Jac*w[j])
	    end # Loop over quadrature points
	    gatherdofnums!(u, dofnums, fes.conn[i]); # retrieve degrees of freedom
	    assemble!(assembler, (elmat + elmat')/2, dofnums, dofnums); # assemble symmetric matrix
	end # Loop over elements
	return makematrix!(assembler);
end

function infsup_gh(self::AbstractFEMMDeforLinear, geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return infsup_gh(self, assembler, geom, u);
end

function _buffers3(self::AbstractFEMMDeforLinear, geom::NodalField, u::NodalField)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare buffers
    elmat = fill(zero(FFlt), elmatdim, elmatdim);      # element matrix -- buffer
    dofnums = zeros(FInt, elmatdim); # degree of freedom array -- buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    gradN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    return dofnums, loc, J, csmatTJ, gradN, elmat
end

"""
    infsup_sh(self::AbstractFEMMDeforLinear, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the matrix to produce the seminorm of the displacement (square root of
the sum of the squares of the derivatives of the components of displacement).

This matrix is used in the numerical infsup test (Klaus-Jurgen Bathe, The
inf-sup condition and its evaluation for mixed finite element methods,
Computers and Structures 79 (2001) 243-252.)

!!! note 
This computation has not been optimized in any way. It can be expected to be
inefficient.
"""
function infsup_sh(self::AbstractFEMMDeforLinear, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
	fes = self.integdomain.fes
	npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
	dofnums, loc, J, csmatTJ, gradN, elmat = _buffers3(self, geom, u)
	startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), u.nfreedofs, u.nfreedofs);
	for i = 1:count(fes) # Loop over elements
	    fill!(elmat,  0.0); # Initialize element matrix
	    for j = 1:npts # Loop over quadrature points
	        locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
	        Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
	        updatecsmat!(self.mcsys, loc, J, fes.label[i]);
	        At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
	        gradN!(fes, gradN, gradNparams[j], csmatTJ);
	        vgradm = vgradmat(self.mr, Ns[j], gradN, loc);
	        elmat += (transpose(vgradm) * vgradm) * (Jac*w[j])
	    end # Loop over quadrature points
	    gatherdofnums!(u, dofnums, fes.conn[i]); # retrieve degrees of freedom
	    assemble!(assembler, (elmat + elmat')/2, dofnums, dofnums); # assemble symmetric matrix
	end # Loop over elements
	return makematrix!(assembler);
end

function infsup_sh(self::AbstractFEMMDeforLinear, geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return infsup_sh(self, assembler, geom, u);
end

end
