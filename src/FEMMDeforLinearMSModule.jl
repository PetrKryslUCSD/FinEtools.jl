"""
    FEMMDeforLinearMSModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear deformation models:  mean-strain  formulation.
"""
module FEMMDeforLinearMSModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: AbstractFESet, FESetH8, FESetT10, manifdim, nodesperelem, gradN!
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FEMMDeforLinearBaseModule: AbstractFEMMDeforLinear
import FinEtools.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D
import FinEtools.MatDeforLinearElasticModule: AbstractMatDeforLinearElastic, tangentmoduli!, update!, thermalstrain!
import FinEtools.MatDeforElastIsoModule: MatDeforElastIso
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!, gathervalues_asmat!
import FinEtools.NodalFieldModule: NodalField
import FinEtools.CSysModule: CSys, updatecsmat!
import FinEtools.DeforModelRedModule: nstressstrain, nthermstrain, Blmat! 
import FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
using FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, add_btv!, loc!, jac!, locjac!
import FinEtools.FEMMDeforLinearBaseModule: stiffness, nzebcloadsstiffness, mass, thermalstrainloads, inspectintegpoints
import FinEtools.FEMMBaseModule: associategeometry!
import FinEtools.MatDeforModule: rotstressvec
import LinearAlgebra: mul!, Transpose, UpperTriangular
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
import LinearAlgebra: norm, qr, diag, dot, cond
import Statistics: mean

"""
    AbstractFEMMDeforLinearMS <: AbstractFEMMDeforLinear

Abstract type for mean-strain linear deformation FEMM.
"""
abstract type AbstractFEMMDeforLinearMS <: AbstractFEMMDeforLinear end

"""
    FEMMDeforLinearMSH8{MR<:AbstractDeforModelRed, S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic} <: AbstractFEMMDeforLinearMS

Type for mean-strain linear deformation FEMM based on eight-node hexahedral elements.
"""
mutable struct FEMMDeforLinearMSH8{MR<:AbstractDeforModelRed, S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic} <: AbstractFEMMDeforLinearMS
    mr::Type{MR}
    integdomain::IntegDomain{S, F} # geometry data
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
    stabilization_material::MatDeforElastIso
    phis::Vector{FFlt}
end

function FEMMDeforLinearMSH8(mr::Type{MR}, integdomain::IntegDomain{S, F}, mcsys::CSys, material::M) where {MR<:AbstractDeforModelRed,  S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearMSH8(mr, integdomain, mcsys, material, stabilization_material, fill(zero(FFlt), 1))
end

function FEMMDeforLinearMSH8(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed,  S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearMSH8(mr, integdomain, CSys(manifdim(integdomain.fes)), material, stabilization_material, fill(zero(FFlt), 1))
end

"""
    FEMMDeforLinearMST10{MR<:AbstractDeforModelRed, S<:FESetT10, F<:Function, M<:AbstractMatDeforLinearElastic} <: AbstractFEMMDeforLinearMS

Type for mean-strain linear deformation FEMM based on 10-node tetrahedral elements.
"""
mutable struct FEMMDeforLinearMST10{MR<:AbstractDeforModelRed, S<:FESetT10, F<:Function, M<:AbstractMatDeforLinearElastic} <: AbstractFEMMDeforLinearMS
    mr::Type{MR}
    integdomain::IntegDomain{S, F} # geometry data
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
    stabilization_material::MatDeforElastIso
    phis::Vector{FFlt}
end

function FEMMDeforLinearMST10(mr::Type{MR}, integdomain::IntegDomain{S, F}, mcsys::CSys, material::M) where {MR<:AbstractDeforModelRed,  S<:FESetT10, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearMST10(mr, integdomain, mcsys, material, stabilization_material, fill(zero(FFlt), 1))
end

function FEMMDeforLinearMST10(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed,  S<:FESetT10, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearMST10(mr, integdomain, CSys(manifdim(integdomain.fes)), material, stabilization_material, fill(zero(FFlt), 1))
end

function _make_stabilization_material(material::M) where {M}
    ns = fieldnames(typeof(material))
    E = 0.0; nu = 0.0
    if :E in ns
        E = material.E
        if material.nu < 0.3
            nu = material.nu
        else
            nu = 0.3 + (material.nu - 0.3) / 2.0
        end
    else
        if :E1 in ns
            E = mean([material.E1, material.E2, material.E3])
            nu = min(material.nu12, material.nu13, material.nu23)
        else
            error("No clues on how to construct the stabilization material")
        end
    end
    return  MatDeforElastIso(material.mr, 0.0, E, nu, 0.0)
end

function _buffers1(self::AbstractFEMMDeforLinearMS, geom::NodalField, npts::FInt)
    fes = self.integdomain.fes
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    # Prepare buffers
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    gradN = fill(zero(FFlt), nne, mdim);
    return loc, J, csmatTJ, gradN
end

function _buffers2(self::AbstractFEMMDeforLinearMS, geom::NodalField, u::NodalField, npts::FInt)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare buffers
    elmat = fill(zero(FFlt), elmatdim, elmatdim);      # element matrix -- buffer
    dofnums = zeros(FInt, 1, elmatdim); # degree of freedom array -- buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    AllgradN = FFltMat[];
    for ixxxx = 1:npts
        push!(AllgradN, fill(zero(FFlt), nne, mdim));
    end
    Jac = fill(zero(FFlt), npts);
    MeangradN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    D = fill(zero(FFlt), nstrs, nstrs); # material stiffness matrix -- buffer
    Dstab = fill(zero(FFlt), nstrs, nstrs); # material stiffness matrix -- buffer
    B = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    DB = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    Bbar = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    elvecfix = fill(zero(FFlt), elmatdim); # vector of prescribed displ. -- buffer
    elvec = fill(zero(FFlt), elmatdim); # element vector -- buffer
    return dofnums, loc, J, csmatTJ, AllgradN, MeangradN, Jac, D, Dstab, B, DB, Bbar, elmat, elvec, elvecfix
end

function centroid!(self::F,  loc, X::FFltMat, conn::C) where {F<:FEMMDeforLinearMSH8, C}
    weights = [0.125
                0.125
                0.125
                0.125
                0.125
                0.125
                0.125
                0.125]
    return loc!(loc, X, conn, reshape(weights, 8, 1))
end

function centroid!(self::F, loc, X::FFltMat, conn::C) where {F<:FEMMDeforLinearMST10, C}
    weights = [ -0.125
                -0.125
                -0.125
                -0.125
                0.250
                0.250
                0.250
                0.250
                0.250
                0.250]
    return loc!(loc, X, conn, reshape(weights, 10, 1)) 
end

"""
    associategeometry!(self::FEMMAbstractBase,  geom::NodalField{FFlt})

Associate geometry field with the FEMM.

Compute the  correction factors to account for  the shape of the  elements.
"""
function associategeometry!(self::F,  geom::NodalField{FFlt}) where {F<:FEMMDeforLinearMSH8}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    loc, J, csmatTJ, gradN = _buffers1(self, geom, npts)
    self.phis = fill(zero(FFlt), count(fes))
    for i = 1:count(fes) # Loop over elements
        # NOTE: the coordinate system should be evaluated at a single point within the
        # element in order for the derivatives to be consistent at all quadrature points
        loc = centroid!(self,  loc, geom.values, fes.conn[i])
        updatecsmat!(self.mcsys, loc, J, fes.label[i]);
        for j = 1:npts # Loop over quadrature points
            jac!(J, geom.values, fes.conn[i], gradNparams[j]) 
            At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], csmatTJ);
            h2 = diag(transpose(csmatTJ)*csmatTJ)
            cap_phi = (2 * (1 + self.stabilization_material.nu) * (minimum(h2) / maximum(h2)))  # Plane stress
            phi = cap_phi / (1 + cap_phi)
            self.phis[i] = max(self.phis[i], phi)
        end # Loop over quadrature points
    end # Loop over elements
    return self
end

"""
    associategeometry!(self::FEMMAbstractBase,  geom::NodalField{FFlt})

Associate geometry field with the FEMM.

Compute the  correction factors to account for  the shape of the  elements.
"""
function associategeometry!(self::F,  geom::NodalField{FFlt}) where {F<:FEMMDeforLinearMST10}
    gamma = 2.6; C = 1e4;
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    loc, J, csmatTJ, gradN = _buffers1(self, geom, npts)
    self.phis = fill(zero(FFlt), count(fes))
    for i = 1:count(fes) # Loop over elements
        # NOTE: the coordinate system should be evaluated at a single point within the
        # element in order for the derivatives to be consistent at all quadrature points
        loc = centroid!(self,  loc, geom.values, fes.conn[i])
        updatecsmat!(self.mcsys, loc, J, fes.label[i]);
        for j = 1:npts # Loop over quadrature points
            jac!(J, geom.values, fes.conn[i], gradNparams[j]) 
            At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            condJ = cond(csmatTJ);
            cap_phi = C*(1.0/condJ)^(gamma);
            phi = cap_phi / (1 + cap_phi)
            self.phis[i] = max(self.phis[i], phi)
        end # Loop over quadrature points
    end # Loop over elements
    return self
end

"""
stiffness(self::AbstractFEMMDeforLinearMS, assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute and assemble  stiffness matrix.
"""
function stiffness(self::AbstractFEMMDeforLinearMS, assembler::A,
    geom::NodalField{FFlt},
    u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, AllgradN, MeangradN, Jac, D, Dstab, B, DB, Bbar, elmat, elvec, elvecfix = _buffers2(self, geom, u, npts)
    realmat = self.material
    stabmat = self.stabilization_material
    tangentmoduli!(realmat, D, 0.0, 0.0, loc, 0)
    tangentmoduli!(stabmat, Dstab, 0.0, 0.0, loc, 0)
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes),
    u.nfreedofs, u.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        # NOTE: the coordinate system should be evaluated at a single point within the
        # element in order for the derivatives to be consistent at all quadrature points
        loc = centroid!(self,  loc, geom.values, fes.conn[i])
        updatecsmat!(self.mcsys, loc, J, fes.label[i]);
        vol = 0.0; # volume of the element
        fill!(MeangradN, 0.0) # mean basis function gradients
        for j = 1:npts # Loop over quadrature points
            jac!(J, geom.values, fes.conn[i], gradNparams[j]) 
            Jac[j] = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, AllgradN[j], gradNparams[j], csmatTJ);
            dvol = Jac[j]*w[j]
            MeangradN .= MeangradN .+ AllgradN[j]*dvol
            vol = vol + dvol
        end # Loop over quadrature points
        MeangradN .= MeangradN/vol
        Blmat!(self.mr, Bbar, Ns[1], MeangradN, loc, self.mcsys.csmat);
        fill!(elmat,  0.0); # Initialize element matrix
        add_btdb_ut_only!(elmat, Bbar, vol, D, DB)
        add_btdb_ut_only!(elmat, Bbar, -self.phis[i]*vol, Dstab, DB)
        for j = 1:npts # Loop over quadrature points
            Blmat!(self.mr, B, Ns[j], AllgradN[j], loc, self.mcsys.csmat);
            add_btdb_ut_only!(elmat, B, self.phis[i]*Jac[j]*w[j], Dstab, DB)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(u, dofnums, fes.conn[i]); # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end


"""
nzebcloadsstiffness(self::AbstractFEMMDeforLinear,  assembler::A,
  geom::NodalField{FFlt},
  u::NodalField{T}) where {A<:AbstractSysvecAssembler, T<:Number}

Compute load vector for nonzero EBC for fixed displacement.
"""
function nzebcloadsstiffness(self::AbstractFEMMDeforLinearMS,  assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysvecAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, AllgradN, MeangradN, Jac, D, Dstab, B, DB, Bbar, elmat, elvec, elvecfix = _buffers2(self, geom, u, npts)
    realmat = self.material
    stabmat = self.stabilization_material
    tangentmoduli!(realmat, D, 0.0, 0.0, loc, 0)
    tangentmoduli!(stabmat, Dstab, 0.0, 0.0, loc, 0)
    startassembly!(assembler,  u.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gatherfixedvalues_asvec!(u, elvecfix, fes.conn[i]);# retrieve element displacement vector
        if norm(elvecfix, Inf) != 0.0   # Is the load nonzero?
            # NOTE: the coordinate system should be evaluated at a single point within the
            # element in order for the derivatives to be consistent at all quadrature points
            loc = centroid!(self,  loc, geom.values, fes.conn[i])
            updatecsmat!(self.mcsys, loc, J, fes.label[i]);
            vol = 0.0; # volume of the element
            fill!(MeangradN, 0.0) # mean basis function gradients
            for j = 1:npts # Loop over quadrature points
                jac!(J, geom.values, fes.conn[i], gradNparams[j]) 
                Jac[j] = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
                At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
                gradN!(fes, AllgradN[j], gradNparams[j], csmatTJ);
                dvol = Jac[j]*w[j]
                MeangradN .= MeangradN .+ AllgradN[j]*dvol
                vol = vol + dvol
            end # Loop over quadrature points
            MeangradN .= MeangradN/vol
            Blmat!(self.mr, Bbar, Ns[1], MeangradN, loc, self.mcsys.csmat);
            fill!(elmat,  0.0); # Initialize element matrix
            add_btdb_ut_only!(elmat, Bbar, vol, D, DB)
            add_btdb_ut_only!(elmat, Bbar, -self.phis[i]*vol, Dstab, DB)
            for j = 1:npts # Loop over quadrature points
                Blmat!(self.mr, B, Ns[j], AllgradN[j], loc, self.mcsys.csmat);
                add_btdb_ut_only!(elmat, B, self.phis[i]*Jac[j]*w[j], Dstab, DB)
            end # Loop over quadrature points
            complete_lt!(elmat)
            gatherdofnums!(u, dofnums, fes.conn[i]); # retrieve degrees of freedom
            A_mul_B!(elvec, elmat, elvecfix)
            assemble!(assembler,  -elvec,  dofnums); # assemble element load vector
        end
    end # Loop over elements
    return makevector!(assembler);
end

function nzebcloadsstiffness(self::AbstractFEMMDeforLinearMS, geom::NodalField{FFlt}, u::NodalField{T}) where {T<:Number}
    assembler = SysvecAssembler()
    return  nzebcloadsstiffness(self, assembler, geom, u);
end

function _iip_meanonly(self::AbstractFEMMDeforLinearMS, geom::NodalField{FFlt},  u::NodalField{T}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F,  idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, AllgradN, MeangradN, Jac, D, Dstab, B, DB, Bbar, elmat, elvec, elvecfix = _buffers2(self, geom, u, npts)
    MeanN = deepcopy(Ns[1])
    realmat = self.material
    stabmat = self.stabilization_material
    # Sort out  the output requirements
    outputcsys = deepcopy(self.mcsys); # default: report the stresses in the material coord system
    for apair in pairs(context)
        sy, val = apair
        if sy == :outputcsys
            outputcsys = deepcopy(val)
        end
    end
    t= 0.0
    dt = 0.0
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    xe = fill(zero(FFlt), nne, sdim); # array of node coordinates -- buffer
    dTe = fill(zero(FFlt), nodesperelem(fes)) # nodal temperatures -- buffer
    ue = fill(zero(FFlt), size(elmat, 1)); # array of node displacements -- buffer
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
        gathervalues_asvec!(dT, dTe, fes.conn[i]);# retrieve element temperature increments
        # NOTE: the coordinate system should be evaluated at a single point within the
        # element in order for the derivatives to be consistent at all quadrature points
        loc = centroid!(self,  loc, geom.values, fes.conn[i])
        updatecsmat!(self.mcsys, loc, J, fes.label[i]);
        updatecsmat!(outputcsys, loc, J, fes.label[i]);
        vol = 0.0; # volume of the element
        fill!(MeangradN, 0.0) # mean basis function gradients
        fill!(MeanN, 0.0) # mean basis function gradients
        for j = 1:npts # Loop over quadrature points
            jac!(J, geom.values, fes.conn[i], gradNparams[j]) 
            Jac[j] = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, AllgradN[j], gradNparams[j], csmatTJ);
            dvol = Jac[j]*w[j]
            MeangradN .= MeangradN .+ AllgradN[j]*dvol
            MeanN .= MeanN .+ Ns[j]*dvol
            vol = vol + dvol
        end # Loop over quadrature points
        MeangradN .= MeangradN/vol
        Blmat!(self.mr, Bbar, MeanN, MeangradN, loc, self.mcsys.csmat);
        MeanN .= MeanN/vol
        qpdT = dot(vec(dTe), vec(MeanN));# Quadrature point temperature increment
        # Quadrature point quantities
        A_mul_B!(qpstrain, Bbar, ue); # strain in material coordinates
        thermalstrain!(realmat, qpthstrain, qpdT)
        # Material updates the state and returns the output
        out = update!(realmat, qpstress, out, vec(qpstrain), qpthstrain, t, dt, loc, fes.label[i], quantity)
        if (quantity == :Cauchy)   # Transform stress tensor,  if that is "quantity"
            (length(out1) >= length(out)) || (out1 = zeros(length(out)))
            rotstressvec(self.mr, out1, out, transpose(self.mcsys.csmat))# To global coord sys
            rotstressvec(self.mr, out, out1, outputcsys.csmat)# To output coord sys
        end
        # Call the inspector
        idat = inspector(idat, i, fes.conn[i], xe, out, loc);
    end # Loop over elements
    return idat; # return the updated inspector data
end

function _iip_extrapmean(self::AbstractFEMMDeforLinearMS, geom::NodalField{FFlt},  u::NodalField{T}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F,  idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, AllgradN, MeangradN, Jac, D, Dstab, B, DB, Bbar, elmat, elvec, elvecfix = _buffers2(self, geom, u, npts)
    MeanN = deepcopy(Ns[1])
    realmat = self.material
    stabmat = self.stabilization_material
    # Sort out  the output requirements
    outputcsys = deepcopy(self.mcsys); # default: report the stresses in the material coord system
    for apair in pairs(context)
        sy, val = apair
        if sy == :outputcsys
            outputcsys = val
        end
    end
    t= 0.0
    dt = 0.0
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    xe = fill(zero(FFlt), nne, sdim); # array of node coordinates -- buffer
    dTe = fill(zero(FFlt), nodesperelem(fes)) # nodal temperatures -- buffer
    ue = fill(zero(FFlt), size(elmat, 1)); # array of node displacements -- buffer
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
        gathervalues_asvec!(dT, dTe, fes.conn[i]);# retrieve element temperature increments
        # NOTE: the coordinate system should be evaluated at a single point within the
        # element in order for the derivatives to be consistent at all quadrature points
        loc = centroid!(self,  loc, geom.values, fes.conn[i])
        updatecsmat!(self.mcsys, loc, J, fes.label[i]);
        updatecsmat!(outputcsys, loc, J, fes.label[i]);
        vol = 0.0; # volume of the element
        fill!(MeangradN, 0.0) # mean basis function gradients
        fill!(MeanN, 0.0) # mean basis function gradients
        for j = 1:npts # Loop over quadrature points
            jac!(J, geom.values, fes.conn[i], gradNparams[j]) 
            Jac[j] = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, AllgradN[j], gradNparams[j], csmatTJ);
            dvol = Jac[j]*w[j]
            MeangradN .= MeangradN .+ AllgradN[j]*dvol
            MeanN .= MeanN .+ Ns[j]*dvol
            vol = vol + dvol
        end # Loop over quadrature points
        MeangradN .= MeangradN/vol
        Blmat!(self.mr, Bbar, MeanN, MeangradN, loc, self.mcsys.csmat);
        MeanN .= MeanN/vol
        qpdT = dot(vec(dTe), vec(MeanN));# Quadrature point temperature increment
        # Quadrature point quantities
        A_mul_B!(qpstrain, Bbar, ue); # strain in material coordinates
        thermalstrain!(realmat, qpthstrain, qpdT)
        # Material updates the state and returns the output
        out = update!(realmat, qpstress, out,
            vec(qpstrain), qpthstrain, t, dt, loc, fes.label[i], quantity)
        if (quantity == :Cauchy)   # Transform stress tensor,  if that is "quantity"
            (length(out1) >= length(out)) || (out1 = zeros(length(out)))
            rotstressvec(self.mr, out1, out, transpose(self.mcsys.csmat))# To global coord sys
            rotstressvec(self.mr, out, out1, outputcsys.csmat)# To output coord sys
        end
        # Call the inspector for each node location
        for nod = 1:size(xe, 1)
            idat = inspector(idat, i, fes.conn[i], xe, out, xe[nod, :]);
        end
    end # Loop over elements
    return idat; # return the updated inspector data
end

function _iip_extraptrend(self::AbstractFEMMDeforLinearMS, geom::NodalField{FFlt},  u::NodalField{T}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F,  idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, AllgradN, MeangradN, Jac, D, Dstab, B, DB, Bbar, elmat, elvec, elvecfix = _buffers2(self, geom, u, npts)
    MeanN = deepcopy(Ns[1])
    realmat = self.material
    stabmat = self.stabilization_material
    # Sort out  the output requirements
    outputcsys = deepcopy(self.mcsys); # default: report the stresses in the material coord system
    for apair in pairs(context)
        sy, val = apair
        if sy == :outputcsys
            outputcsys = val
        end
    end
    t= 0.0
    dt = 0.0
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    xe = fill(zero(FFlt), nne, sdim); # array of node coordinates -- buffer
    dTe = fill(zero(FFlt), nodesperelem(fes)) # nodal temperatures -- buffer
    ue = fill(zero(FFlt), size(elmat, 1)); # array of node displacements -- buffer
    qpdT = 0.0; # node temperature increment
    qpstrain = fill(zero(FFlt), nstressstrain(self.mr), 1); # total strain -- buffer
    qpthstrain = fill(zero(FFlt), nthermstrain(self.mr)); # thermal strain -- buffer
    qpstress = fill(zero(FFlt), nstressstrain(self.mr)); # stress -- buffer
    rout1 = fill(zero(FFlt), nstressstrain(self.mr)); # stress -- buffer
    rout =  fill(zero(FFlt), nstressstrain(self.mr));# output -- buffer
    sbout = deepcopy(rout)
    sbout1 = deepcopy(sbout)
    sout = deepcopy(rout)
    sout1 = deepcopy(sout)
    sqploc = deepcopy(loc)
    sstoredout = fill(zero(FFlt), npts, length(sout))
    A = ones(FFlt, npts, 4)
    # Loop over  all the elements and all the quadrature points within them
    for ilist = 1:length(felist) # Loop over elements
        i = felist[ilist];
        gathervalues_asmat!(geom, xe, fes.conn[i]);# retrieve element coords
        gathervalues_asvec!(u, ue, fes.conn[i]);# retrieve element displacements
        gathervalues_asvec!(dT, dTe, fes.conn[i]);# retrieve element temperature increments
        # NOTE: the coordinate system should be evaluated at a single point within the
        # element in order for the derivatives to be consistent at all quadrature points
        loc = centroid!(self,  loc, geom.values, fes.conn[i])
        updatecsmat!(self.mcsys, loc, J, fes.label[i]);
        updatecsmat!(outputcsys, loc, J, fes.label[i]);
        vol = 0.0; # volume of the element
        fill!(MeangradN, 0.0) # mean basis function gradients
        fill!(MeanN, 0.0) # mean basis function gradients
        for j = 1:npts # Loop over quadrature points
            jac!(J, geom.values, fes.conn[i], gradNparams[j]) 
            Jac[j] = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, AllgradN[j], gradNparams[j], csmatTJ);
            dvol = Jac[j]*w[j]
            MeangradN .= MeangradN .+ AllgradN[j]*dvol
            MeanN .= MeanN .+ Ns[j]*dvol
            vol = vol + dvol
        end # Loop over quadrature points
        MeangradN .= MeangradN/vol
        Blmat!(self.mr, Bbar, MeanN, MeangradN, loc, self.mcsys.csmat);
        MeanN .= MeanN/vol
        qpdT = dot(vec(dTe), vec(MeanN));# Quadrature point temperature increment
        # Quadrature point quantities
        A_mul_B!(qpstrain, Bbar, ue); # strain in material coordinates
        thermalstrain!(realmat, qpthstrain, qpdT)
        # REAL Material updates the state and returns the output
        rout = update!(realmat, qpstress, rout, vec(qpstrain), qpthstrain, t, dt, loc, fes.label[i], quantity)
        if (quantity == :Cauchy)   # Transform stress tensor,  if that is "quantity"
            (length(rout1) >= length(rout)) || (rout1 = zeros(length(rout)))
            rotstressvec(self.mr, rout1, rout, transpose(self.mcsys.csmat))# To global coord sys
            rotstressvec(self.mr, rout, rout1, outputcsys.csmat)# To output coord sys
        end
        # STABILIZATION Material updates the state and returns the output
        sbout = update!(stabmat, qpstress, sbout, vec(qpstrain), qpthstrain, t, dt, loc, fes.label[i], quantity)
        if (quantity == :Cauchy)   # Transform stress tensor,  if that is "quantity"
            (length(sbout1) >= length(sbout)) || (sbout1 = zeros(length(sbout)))
            rotstressvec(self.mr, sbout1, sbout, transpose(self.mcsys.csmat))# To global coord sys
            rotstressvec(self.mr, sbout, sbout1, outputcsys.csmat)# To output coord sys
        end
        for j = 1:npts # Loop over quadrature points (STABILIZATION material)
            At_mul_B!(sqploc, Ns[j], xe);# Quadrature point location
            A[j, 1:3] .= vec(sqploc - loc);
            Blmat!(self.mr, B, Ns[j], AllgradN[j], sqploc, self.mcsys.csmat);
            qpdT = dot(vec(dTe), vec(Ns[j]));# Quadrature point temperature increment
            #  Quadrature point quantities
            A_mul_B!(qpstrain, B, ue); # strain in material coordinates
            thermalstrain!(stabmat, qpthstrain, qpdT)
            # Material updates the state and returns the output
            sout = update!(stabmat, qpstress, sout, vec(qpstrain), qpthstrain, t, dt, loc, fes.label[i], quantity)
            if (quantity == :Cauchy)   # Transform stress tensor,  if that is "quantity"
                (length(sout1) >= length(sout)) || (sout1 = zeros(length(sout)))
                rotstressvec(self.mr, sout1, sout, transpose(self.mcsys.csmat))# To global coord sys
                rotstressvec(self.mr, sout, sout1, outputcsys.csmat)# To output coord sys
            end
            sstoredout[j, :] .= sout # store  the output for this quadrature point
        end # Loop over quadrature points
        #  Solve for the least-square fit parameters
        fact = qr(A); 
        p = UpperTriangular(fact.R)\(transpose(Array(fact.Q))*sstoredout) # R \ (transpose(Q) * sstoredout)
        for nod = 1:size(xe, 1)
            #  Predict the value  of the output quantity at the node
            xdel = vec(@view xe[nod, :]) - vec(loc)
            nout = rout + self.phis[i] * (- sbout + vec(reshape(xdel, 1, 3) * p[1:3, :]) + p[4, :])
            # Call the inspector for the node location
            idat = inspector(idat, i, fes.conn[i], xe, nout, xe[nod, :]);
        end
    end # Loop over elements
    return idat; # return the updated inspector data
end

"""
   inspectintegpoints(self::FEMMDeforLinear,
     geom::NodalField{FFlt},  u::NodalField{T},
     dT::NodalField{FFlt},
     felist::FIntVec,
     inspector::F,  idat, quantity=:Cauchy;
     context...) where {T<:Number, F<:Function}

Inspect integration point quantities.

`geom` - reference geometry field
`u` - displacement field
`dT` - temperature difference field
`felist` - indexes of the finite elements that are to be inspected:
    The fes to be included are: `fes[felist]`.
`context`    - structure: see the update!() method of the material.
`inspector` - functionwith the signature
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
function inspectintegpoints(self::AbstractFEMMDeforLinearMS, geom::NodalField{FFlt},  u::NodalField{T}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F,  idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    reportat = :meanonly
    for apair in pairs(context)
        sy, val = apair
        if sy == :reportat
            reportat = val
        end
    end
    if reportat == :extraptrend
        return _iip_extraptrend(self, geom, u, dT, felist,
            inspector, idat, quantity; context...);
    elseif reportat == :extrapmean
        return _iip_extrapmean(self, geom, u, dT, felist,
            inspector, idat, quantity; context...);
    elseif reportat == :meanonly || true # this is the default
        return _iip_meanonly(self, geom, u, dT, felist,
            inspector, idat, quantity; context...);
    end
    return idat
end

function inspectintegpoints(self::AbstractFEMMDeforLinearMS, geom::NodalField{FFlt},  u::NodalField{T}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    dT = NodalField(fill(zero(FFlt), nnodes(geom), 1)) # zero difference in temperature
    return inspectintegpoints(self, geom, u, dT, felist, inspector, idat, quantity; context...);
end

end
