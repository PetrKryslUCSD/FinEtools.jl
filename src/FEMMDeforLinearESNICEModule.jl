"""
    FEMMDeforLinearESNICEModule

Formulation for the small displacement, small strain deformation
model for Nodally-Integrated Continuum Elements (NICE).

The approximation is  originally from Dohrmann et al IJNME 47 (2000).
The formulation was subsequently developed in Krysl, P. and Zhu, B.
Locking-free continuum displacement finite elements with nodal
integration, International Journal for Numerical Methods in Engineering,
76,7,1020-1043,2008.

The stabilization scheme comes from papers on energy-sampling stabilization
for mean-strain elements (Krysl).
"""
module FEMMDeforLinearESNICEModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: AbstractFESet, FESetH8, FESetT4, manifdim, nodesperelem, gradN!
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FEMMDeforLinearBaseModule: AbstractFEMMDeforLinear
import FinEtools.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D
import FinEtools.MatDeforLinearElasticModule: AbstractMatDeforLinearElastic, tangentmoduli!, update!, thermalstrain!
import FinEtools.MatDeforElastIsoModule: MatDeforElastIso
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!, gathervalues_asmat!
import FinEtools.NodalFieldModule: NodalField, nnodes
import FinEtools.CSysModule: CSys, updatecsmat!
import FinEtools.FENodeToFEMapModule: FENodeToFEMap
import FinEtools.DeforModelRedModule: nstressstrain, nthermstrain, Blmat!, divmat, vgradmat
import FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
using FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, add_btv!, loc!, jac!, locjac!, adjugate3!
import FinEtools.FEMMDeforLinearBaseModule: stiffness, nzebcloadsstiffness, mass, thermalstrainloads, inspectintegpoints
import FinEtools.FEMMBaseModule: associategeometry!
import FinEtools.MatDeforModule: rotstressvec!
import FinEtools.MatModule: massdensity
import LinearAlgebra: mul!, Transpose, UpperTriangular, eigvals, det
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
import LinearAlgebra: norm, qr, diag, dot, cond, I, cross
import Statistics: mean
import StatsBase: geomean

"""
    AbstractFEMMDeforLinearESNICE <: AbstractFEMMDeforLinear

Abstract FEMM type for Nodally Integrated Continuum Elements (NICE) with energy-sampling stabilization.
"""
abstract type AbstractFEMMDeforLinearESNICE <: AbstractFEMMDeforLinear end

# Tetrahedron
# The coefficient set below was obtained by fitting the ratio of energies
# true/approximate for the finite element model of six tetrahedra arranged
# into a rectangular block and subject to pure bending
# Fitting for a small aspect-ratio range (1.0 to 10)
_T4_stabilization_parameters = (2.101588423297799,  1.311321055432958)


mutable struct _NodalBasisFunctionGradients
    gradN::FFltMat
    patchconn::FIntVec
    Vpatch::FFlt
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

"""
    FEMMDeforLinearESNICET4{MR<:AbstractDeforModelRed, S<:FESetT4, F<:Function, M<:AbstractMatDeforLinearElastic} <: AbstractFEMMDeforLinearESNICE

FEMM type for Nodally Integrated Continuum Elements (NICE) based on the 4-node tetrahedron.
"""
mutable struct FEMMDeforLinearESNICET4{MR<:AbstractDeforModelRed, S<:FESetT4, F<:Function, M<:AbstractMatDeforLinearElastic} <: AbstractFEMMDeforLinearESNICE
    mr::Type{MR}
    integdomain::IntegDomain{S, F} # geometry data
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
    stabilization_material::MatDeforElastIso
    nodalbasisfunctiongrad::Vector{_NodalBasisFunctionGradients}
    ephis::Vector{FFlt}
    nphis::Vector{FFlt}
end

"""
    FEMMDeforLinearESNICEH8{MR<:AbstractDeforModelRed, S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic} <: AbstractFEMMDeforLinearESNICE

FEMM type for Nodally Integrated Continuum Elements (NICE) based on the 8-node hexahedron.
"""
mutable struct FEMMDeforLinearESNICEH8{MR<:AbstractDeforModelRed, S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic} <: AbstractFEMMDeforLinearESNICE
    mr::Type{MR}
    integdomain::IntegDomain{S, F} # geometry data
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
    stabilization_material::MatDeforElastIso
    nodalbasisfunctiongrad::Vector{_NodalBasisFunctionGradients}
    ephis::Vector{FFlt}
    nphis::Vector{FFlt}
end

function FEMMDeforLinearESNICET4(mr::Type{MR}, integdomain::IntegDomain{S, F}, mcsys::CSys, material::M) where {MR<:AbstractDeforModelRed,  S<:FESetT4, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearESNICET4(mr, integdomain, mcsys, material, stabilization_material, _NodalBasisFunctionGradients[], fill(zero(FFlt), 1), fill(zero(FFlt), 1))
end

function FEMMDeforLinearESNICET4(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed,  S<:FESetT4, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearESNICET4(mr, integdomain, CSys(manifdim(integdomain.fes)), material, stabilization_material, _NodalBasisFunctionGradients[], fill(zero(FFlt), 1), fill(zero(FFlt), 1))
end

function FEMMDeforLinearESNICEH8(mr::Type{MR}, integdomain::IntegDomain{S, F}, mcsys::CSys, material::M) where {MR<:AbstractDeforModelRed,  S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearESNICEH8(mr, integdomain, mcsys, material, stabilization_material, _NodalBasisFunctionGradients[], fill(zero(FFlt), 1), fill(zero(FFlt), 1))
end

function FEMMDeforLinearESNICEH8(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed,  S<:FESetH8, F<:Function, M<:AbstractMatDeforLinearElastic}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (mr == DeforModelRed3D) "3D model required"
    stabilization_material = _make_stabilization_material(material)
    return FEMMDeforLinearESNICEH8(mr, integdomain, CSys(manifdim(integdomain.fes)), material, stabilization_material, _NodalBasisFunctionGradients[], fill(zero(FFlt), 1), fill(zero(FFlt), 1))
end

function _buffers1(self::AbstractFEMMDeforLinearESNICE, geom::NodalField)
    fes = self.integdomain.fes
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    # Prepare buffers
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    adjJ = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    gradN = fill(zero(FFlt), nne, mdim);
    xl = fill(zero(FFlt), nne, mdim);
    lconn = collect(1:nne)
    return loc, J, adjJ, csmatTJ, gradN, xl, lconn
end

function _buffers2(self::AbstractFEMMDeforLinearESNICE, geom::NodalField, u::NodalField, npts::FInt)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare buffers
    elmat = fill(zero(FFlt), elmatdim, elmatdim);      # element matrix -- buffer
    B = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    Jac = fill(zero(FFlt), npts);
    D = fill(zero(FFlt), nstrs, nstrs); # material stiffness matrix -- buffer
    Dstab = fill(zero(FFlt), nstrs, nstrs); # material stiffness matrix -- buffer
    return dofnums, loc, J, csmatTJ, Jac, D, Dstab, elmat, B
end

function _buffers3(self::AbstractFEMMDeforLinearESNICE, geom::NodalField, u::NodalField)
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
    B = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    DB = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    elvecfix = fill(zero(FFlt), elmatdim); # vector of prescribed displ. -- buffer
    elvec = fill(zero(FFlt), elmatdim); # element vector -- buffer
    gradN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    return dofnums, B, DB, elmat, elvec, elvecfix, gradN
end

function _patchconn(fes, gl, thisnn)
    # Generate patch connectivity for a given node (thisnn)
    # from the connectivities of the finite elements attached to it.
    return vcat(collect(setdiff(Set([i for j=1:length(gl) for i in fes.conn[gl[j]]]), thisnn)), [thisnn])
end

function _computenodalbfungrads(self, geom)
    # # Compute the nodal basis function gradients.
    # # Return the cell array of structures with attributes
    # %      bfun_gradients{nix}.Nspd= basis function gradient matrix
    # #        bfun_gradients{nix}.Vpatch= nodal patch volume
    # #        bfun_gradients{nix}.patchconn= nodal patch connectivity

    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    loc, J, adjJ, csmatTJ, gradN, xl, lconn = _buffers1(self, geom)

    # Get the inverse map from finite element nodes to geometric cells
    fen2fe = FENodeToFEMap(fes.conn, nnodes(geom));
    # Initialize the nodal gradients, nodal patch, and patch connectivity
    bfungrads = fill(_NodalBasisFunctionGradients(fill(0.0, 0, 0), fill(0, 0), 0.0), nnodes(geom));
    # Now loop over all finite element nodes in the map
    lnmap = fill(0, length(fen2fe.map)); # Local node map: buffer to speed up operations
    for nix = 1:length(fen2fe.map)
        gl = fen2fe.map[nix];
        thisnn = nix; # We are at this node
        if !isempty(gl) # This node has an element patch in this block
            # establish local numbering of all nodes of the patch @ node thisnn
            p = _patchconn(fes, gl, thisnn);
            np = length(p);
            lnmap[p] .= 1:np;# now store the local numbers
            c = reshape(geom.values[thisnn, :], 1, ndofs(geom))
            updatecsmat!(self.mcsys, c, J, 0);
            gradNavg = fill(0.0, np, ndofs(geom));# preallocate strain-displacement matrix
            Vpatch = 0.0;
            for k = 1:length(gl)
                i = gl[k]
                kconn = collect(fes.conn[i]);
                pci = findfirst(cx -> cx == thisnn, kconn);# at which node in the element are we with this quadrature point?
                @assert 1 <= pci <= nodesperelem(fes)
                # centered coordinates of nodes in the material coordinate system
                for cn = 1:length(kconn)
                    xl[cn, :] = (reshape(geom.values[kconn[cn], :], 1, ndofs(geom)) - c) * self.mcsys.csmat
                end
                jac!(J, xl, lconn, gradNparams[pci])
                At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
                Jac = Jacobianvolume(self.integdomain, J, c, fes.conn[i], Ns[pci]);
                Vpatch += Jac * w[pci];
                sgradN = gradNparams[pci] * adjugate3!(adjJ, J);
                gradNavg[lnmap[kconn],:] += (w[pci] .* sgradN);
            end
            @assert Vpatch != 0
            gradNavg ./= Vpatch;
            bfungrads[nix] = _NodalBasisFunctionGradients(gradNavg, p, Vpatch);
            lnmap[p] .= 0; # Restore the buffer to pristine condition
        end
    end
    self.nodalbasisfunctiongrad = bfungrads
    return self
end

function _tetaspectratiovol(X)
    edge1 = vec(X[2, :] - X[1, :])
    edge2 = vec(X[3, :] - X[1, :])
    edge3 = vec(X[4, :] - X[1, :])
    edge4 = vec(X[3, :] - X[2, :])
    edge5 = vec(X[4, :] - X[3, :])
    edge6 = vec(X[4, :] - X[2, :])
    V = dot(edge3, cross(edge1, edge2))
    A1 = norm(cross(edge1, edge2)) # This is twice the area of the triangle
    A2 = norm(cross(edge2, edge3))
    A3 = norm(cross(edge3, edge1))
    A4 = norm(cross(edge4, edge6))
    h1, h2, h3, h4 = V/A1, V/A2, V/A3, V/A4
    L1, L2, L3, L4, L5, L6 = norm(edge1), norm(edge2), norm(edge3), norm(edge4), norm(edge5), norm(edge6)
    f = maximum
    return h1/f([L1, L2, L4]), h2/f([L3, L2, L5]), h3/f([L1, L3, L6]), h4/f([L6, L5, L4]), V/6
end

"""
    associategeometry!(self::F,  geom::NodalField{FFlt}) where {F<:FEMMDeforLinearESNICET4}

Associate geometry field with the FEMM.

Compute the  correction factors to account for  the shape of the  elements.
"""
function associategeometry!(self::F,  geom::NodalField{FFlt}; stabilization_parameters = _T4_stabilization_parameters) where {F<:FEMMDeforLinearESNICET4}
    (a, b) = stabilization_parameters
    fes = self.integdomain.fes
    self.ephis = fill(zero(FFlt), count(fes))
    evols = fill(zero(FFlt), count(fes))
    self.nphis = fill(zero(FFlt), nnodes(geom))
    nvols = fill(zero(FFlt), nnodes(geom))
    for i = 1:count(fes) # Loop over elements
        ar1, ar2, ar3, ar4, V = _tetaspectratiovol(geom.values[collect(fes.conn[i]), :])
        evols[i] = V;
        self.ephis[i] = (1.0 / (b * min(ar1, ar2, ar3, ar4) ^a) + 1.0) ^(-1)
        # Accumulate: the stabilization factor at the node is the weighted mean of the stabilization factors of the elements at that node
        for k = 1:nodesperelem(fes)
            nvols[fes.conn[i][k]] += evols[i]
            self.nphis[fes.conn[i][k]] += self.ephis[i] * evols[i]
        end
    end # Loop over elements
    # Now scale the values at the nodes with the nodal volumes
    for k = 1:length(nvols)
        self.nphis[k] /= nvols[k]
    end
    # Now calculate the nodal basis function gradients
    return _computenodalbfungrads(self, geom)
end

"""
    associategeometry!(self::F,  geom::NodalField{FFlt}) where {F<:FEMMDeforLinearESNICEH8}

Associate geometry field with the FEMM.

Compute the  correction factors to account for  the shape of the  elements.
"""
function associategeometry!(self::F,  geom::NodalField{FFlt}) where {F<:FEMMDeforLinearESNICEH8}
    fes = self.integdomain.fes
    self.ephis = fill(zero(FFlt), count(fes))
    evols = fill(zero(FFlt), count(fes))
    self.nphis = fill(zero(FFlt), nnodes(geom))
    nvols = fill(zero(FFlt), nnodes(geom))
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain);

    for i = 1:count(fes) # Loop over elements
        X = geom.values[collect(fes.conn[i]), :]
        V = 0.0;
        for j = 1:npts
            J = X' * gradNparams[j]
            Jac = det(J)
            @assert Jac > 0 "Nonpositive Jacobian"
            V = V + Jac * w[j]
            hs = sum(J .* J, dims=1)
            phi = 2*(1.0 + self.stabilization_material.nu)*minimum(hs)/maximum(hs)
            self.ephis[i] = max(self.ephis[i], phi/(1+phi))
        end
        evols[i] = V;
        # Accumulate: the stabilization factor at the node is the weighted mean of the stabilization factors of the elements at that node
        for k = 1:nodesperelem(fes)
            nvols[fes.conn[i][k]] += evols[i]
            self.nphis[fes.conn[i][k]] += self.ephis[i] * evols[i]
        end
    end # Loop over elements
    # Now scale the values at the nodes with the nodal volumes
    for k = 1:length(nvols)
        self.nphis[k] /= nvols[k]
    end
    # Now calculate the nodal basis function gradients
    return _computenodalbfungrads(self, geom)
end
"""
    stiffness(self::AbstractFEMMDeforLinearESNICE, assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute and assemble  stiffness matrix.
"""
function stiffness(self::AbstractFEMMDeforLinearESNICE, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, Jac, D, Dstab = _buffers2(self, geom, u, npts)
    tangentmoduli!(self.material, D, 0.0, 0.0, loc, 0)
    tangentmoduli!(self.stabilization_material, Dstab, 0.0, 0.0, loc, 0)
    elmatsizeguess = 4*nodesperelem(fes)*ndofs(u)
    startassembly!(assembler, elmatsizeguess, elmatsizeguess, nnodes(u) + count(fes), u.nfreedofs, u.nfreedofs);
    for nix = 1:length(self.nodalbasisfunctiongrad)
        gradN = self.nodalbasisfunctiongrad[nix].gradN
        patchconn = self.nodalbasisfunctiongrad[nix].patchconn
        Vpatch = self.nodalbasisfunctiongrad[nix].Vpatch
        c = reshape(geom.values[nix, :], 1, ndofs(geom))
        updatecsmat!(self.mcsys, c, J, 0);
        nd = length(patchconn) * ndofs(u)
        Bnodal = fill(0.0, size(D, 1), nd)
        Blmat!(self.mr, Bnodal, Ns[1], gradN, c, self.mcsys.csmat);
        elmat = fill(0.0, nd, nd) # Can we SPEED it UP?
        DB = fill(0.0, size(D, 1), nd)
        add_btdb_ut_only!(elmat, Bnodal, Vpatch, D, DB)
        add_btdb_ut_only!(elmat, Bnodal, -self.nphis[nix]*Vpatch, Dstab, DB)
        complete_lt!(elmat)
        dofnums = fill(0, nd)
        gatherdofnums!(u, dofnums, patchconn); # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); # assemble symmetric matrix
    end # Loop over elements
    dofnums, B, DB, elmat, elvec, elvecfix, gradN = _buffers3(self, geom, u)
    # OPTIMIZATION: switch to a single-point quadrature rule here
    for i = 1:count(fes) # Loop over elements
        fill!(elmat,  0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]);
            At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], csmatTJ);
            Blmat!(self.mr, B, Ns[j], gradN, loc, self.mcsys.csmat);
            add_btdb_ut_only!(elmat, B, self.ephis[i]*Jac*w[j], Dstab, DB)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(u, dofnums, fes.conn[i]); # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function stiffness(self::AbstractFEMMDeforLinearESNICE, geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom, u);
end


"""
nzebcloadsstiffness(self::AbstractFEMMDeforLinear,  assembler::A,
  geom::NodalField{FFlt},
  u::NodalField{T}) where {A<:AbstractSysvecAssembler, T<:Number}

Compute load vector for nonzero EBC for fixed displacement.
"""
function nzebcloadsstiffness(self::AbstractFEMMDeforLinearESNICE,  assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysvecAssembler, T<:Number}
    error("Not implemented yet")
end

function nzebcloadsstiffness(self::AbstractFEMMDeforLinearESNICE, geom::NodalField{FFlt}, u::NodalField{T}) where {T<:Number}
    assembler = SysvecAssembler()
    return  nzebcloadsstiffness(self, assembler, geom, u);
end

"""
   inspectintegpoints(self::AbstractFEMMDeforLinearESNICE, geom::NodalField{FFlt},  u::NodalField{T}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F,  idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}

Inspect integration point quantities.

# Arguments
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
          
# Return
The updated inspector data is returned.
"""
function inspectintegpoints(self::AbstractFEMMDeforLinearESNICE, geom::NodalField{FFlt},  u::NodalField{T}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F,  idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, Jac, D, Dstab = _buffers2(self, geom, u, npts)
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
    qpdT = 0.0; # node temperature increment
    qpstrain = fill(zero(FFlt), nstressstrain(self.mr), 1); # total strain -- buffer
    qpthstrain = fill(zero(FFlt), nthermstrain(self.mr)); # thermal strain -- buffer
    qpstress = fill(zero(FFlt), nstressstrain(self.mr)); # stress -- buffer
    out1 = fill(zero(FFlt), nstressstrain(self.mr)); # stress -- buffer
    out =  fill(zero(FFlt), nstressstrain(self.mr));# output -- buffer
    outtot =  fill(zero(FFlt), nstressstrain(self.mr));# output -- buffer
    xe = fill(0.0, nodesperelem(fes), ndofs(geom))
    eue = fill(zero(T), nodesperelem(fes) * ndofs(u))
    gradN = fill(zero(FFlt), nodesperelem(fes), ndofs(u)); 
    B = fill(zero(FFlt), nstressstrain(self.mr), nodesperelem(fes) * ndofs(u)); 
    # Loop over  all the elements and all the quadrature points within them
    for ilist = 1:length(felist) # Loop over elements
    	i = felist[ilist];
    	gathervalues_asmat!(geom, xe, fes.conn[i]);# retrieve element coords
    	for nix = fes.conn[i] # For all nodes connected by this element
    		nodalgradN = self.nodalbasisfunctiongrad[nix].gradN
    		patchconn = self.nodalbasisfunctiongrad[nix].patchconn
    		Vpatch = self.nodalbasisfunctiongrad[nix].Vpatch
    		ue = fill(zero(T), length(patchconn) * ndofs(u))
    		gathervalues_asvec!(u, ue, patchconn);# retrieve element displacements
    		loc = reshape(geom.values[nix, :], 1, ndofs(geom))
    		updatecsmat!(self.mcsys, loc, J, nix);
    		nd = length(patchconn) * ndofs(u)
    		Bnodal = fill(0.0, size(D, 1), nd)
    		Blmat!(self.mr, Bnodal, Ns[1], nodalgradN, loc, self.mcsys.csmat);
    		updatecsmat!(outputcsys, loc, J, nix); # Update output coordinate system
    		# Quadrature point quantities
    		A_mul_B!(qpstrain, Bnodal, ue); # strain in material coordinates
    		qpdT = dT.values[nix] # Quadrature point temperature increment
    		thermalstrain!(self.material, qpthstrain, qpdT)
    		# Real Material: update the state and return the output
    		out = update!(self.material, qpstress, out, vec(qpstrain), qpthstrain, t, dt, loc, nix, quantity)
    		copyto!(outtot, out)
    		# Stabilization Material: update the state and return the output
    		out = update!(self.stabilization_material, qpstress, out, vec(qpstrain), qpthstrain, t, dt, loc, nix, quantity)
    		outtot .+= -self.nphis[nix].*out
    		pci = findfirst(cx -> cx == nix, fes.conn[i]);# at which node are we?
    		locjac!(loc, J, geom.values, fes.conn[i], Ns[pci], gradNparams[pci])
    		updatecsmat!(self.mcsys, loc, J, fes.label[i]);
    		At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
    		gradN!(fes, gradN, gradNparams[pci], csmatTJ);
    		Blmat!(self.mr, B, Ns[pci], gradN, loc, self.mcsys.csmat);
    		gathervalues_asvec!(u, eue, fes.conn[i]);# retrieve element displacements
    		A_mul_B!(qpstrain, B, eue); # strain in material coordinates
    		out = update!(self.stabilization_material, qpstress, out, vec(qpstrain), qpthstrain, t, dt, loc, nix, quantity)
    		outtot .+= +self.ephis[i].*out
    		out, outtot = outtot, out # swap in the total output
    		if (quantity == :Cauchy)   # Transform stress tensor,  if that is "out"
    		    (length(out1) >= length(out)) || (out1 = zeros(length(out)))
    		    rotstressvec!(self.mr, out1, out, transpose(self.mcsys.csmat))# To global coord sys
    		    rotstressvec!(self.mr, out, out1, outputcsys.csmat)# To output coord sys
    		end
    		# Call the inspector
    		idat = inspector(idat, i, fes.conn[i], xe, out, loc);
    	end # Loop over nodes
    end # Loop over elements
    return idat
end

function inspectintegpoints(self::AbstractFEMMDeforLinearESNICE, geom::NodalField{FFlt},  u::NodalField{T}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {T<:Number, F<:Function}
    dT = NodalField(fill(zero(FFlt), nnodes(geom), 1)) # zero difference in temperature
    return inspectintegpoints(self, geom, u, dT, felist, inspector, idat, quantity; context...);
end

"""
    infsup_gh(self::AbstractFEMMDeforLinearESNICE, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the matrix to produce the norm of the divergence of the displacement.

This matrix is used in the numerical infsup test (Klaus-Jurgen Bathe, The
inf-sup condition and its evaluation for mixed finite element methods,
Computers and Structures 79 (2001) 243-252.)

!!! note 
This computation has not been optimized in any way. It can be expected to be
inefficient.
"""
function infsup_gh(self::AbstractFEMMDeforLinearESNICE, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
	fes = self.integdomain.fes
	npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
	elmatsizeguess = 4*nodesperelem(fes)*ndofs(u)
	startassembly!(assembler, elmatsizeguess, elmatsizeguess, nnodes(u) + count(fes), u.nfreedofs, u.nfreedofs);
	for nix = 1:length(self.nodalbasisfunctiongrad)
	    gradN = self.nodalbasisfunctiongrad[nix].gradN
	    patchconn = self.nodalbasisfunctiongrad[nix].patchconn
	    Vpatch = self.nodalbasisfunctiongrad[nix].Vpatch
	    c = reshape(geom.values[nix, :], 1, ndofs(geom))
	    nd = length(patchconn) * ndofs(u)
	    divm = divmat(self.mr, Ns[1], gradN, c);
	    elmat = (transpose(divm) * divm) * Vpatch
	    dofnums = fill(0, nd)
	    gatherdofnums!(u, dofnums, patchconn); # retrieve degrees of freedom
	    assemble!(assembler, elmat, dofnums, dofnums); # assemble symmetric matrix
	end # Loop over nodes
	return makematrix!(assembler);
end

function infsup_gh(self::AbstractFEMMDeforLinearESNICE, geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return infsup_gh(self, assembler, geom, u);
end

"""
    infsup_sh(self::AbstractFEMMDeforLinearESNICE, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the matrix to produce the seminorm of the displacement (square root of
the sum of the squares of the derivatives of the components of displacement).

This matrix is used in the numerical infsup test (Klaus-Jurgen Bathe, The
inf-sup condition and its evaluation for mixed finite element methods,
Computers and Structures 79 (2001) 243-252.)

!!! note 
This computation has not been optimized in any way. It can be expected to be
inefficient.
"""
function infsup_sh(self::AbstractFEMMDeforLinearESNICE, assembler::A, geom::NodalField{FFlt}, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}
	fes = self.integdomain.fes
	npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
	elmatsizeguess = 4*nodesperelem(fes)*ndofs(u)
	startassembly!(assembler, elmatsizeguess, elmatsizeguess, nnodes(u) + count(fes), u.nfreedofs, u.nfreedofs);
	for nix = 1:length(self.nodalbasisfunctiongrad)
	    gradN = self.nodalbasisfunctiongrad[nix].gradN
	    patchconn = self.nodalbasisfunctiongrad[nix].patchconn
	    Vpatch = self.nodalbasisfunctiongrad[nix].Vpatch
	    c = reshape(geom.values[nix, :], 1, ndofs(geom))
	    nd = length(patchconn) * ndofs(u)
	    vgradm = vgradmat(self.mr, Ns[1], gradN, c);
	    elmat = (transpose(vgradm) * vgradm) * Vpatch
	    dofnums = fill(0, nd)
	    gatherdofnums!(u, dofnums, patchconn); # retrieve degrees of freedom
	    assemble!(assembler, elmat, dofnums, dofnums); # assemble symmetric matrix
	end # Loop over nodes
	return makematrix!(assembler);
end

function infsup_sh(self::AbstractFEMMDeforLinearESNICE, geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return infsup_sh(self, assembler, geom, u);
end

end
