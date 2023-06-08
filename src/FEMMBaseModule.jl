"""
    FEMMBaseModule

Module for comments/base operations on interiors and boundaries of domains.
"""
module FEMMBaseModule

__precompile__(true)

import LinearAlgebra: mul!, Transpose
my_At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
import SparseArrays: sparse, findnz
import LinearAlgebra: norm
import ..FENodeSetModule: FENodeSet
import ..FESetModule:
    AbstractFESet,
    manifdim,
    nodesperelem,
    subset,
    map2parametric,
    inparametric,
    centroidparametric,
    bfun
import ..IntegDomainModule: IntegDomain, integrationdata, Jacobianmdim, Jacobianvolume
import ..CSysModule: CSys
import ..FieldModule: ndofs, nents, gatherdofnums!, gathervalues_asmat!
import ..NodalFieldModule: NodalField, nnodes
import ..ElementalFieldModule: ElementalField, nelems
import ..ForceIntensityModule: ForceIntensity, updateforce!
import ..MatrixUtilityModule: locjac!
import ..BoxModule: initbox!, boundingbox, inflatebox!
import ..MeshModificationModule: nodepartitioning, compactnodes, renumberconn!
import ..MeshSelectionModule: selectelem, vselect, findunconnnodes, connectednodes
import ..AssemblyModule:
    AbstractSysvecAssembler,
    AbstractSysmatAssembler,
    SysmatAssemblerSparseSymm,
    startassembly!,
    assemble!,
    makematrix!,
    makevector!,
    SysvecAssembler
import ..FENodeToFEMapModule: FENodeToFEMap

"""
    AbstractFEMM

Abstract type for all finite element model machines.
"""
abstract type AbstractFEMM end

"""
    FEMMBase{S<:AbstractFESet, F<:Function} <: AbstractFEMM

Class for base finite element modeling machine.
"""
mutable struct FEMMBase{ID<:IntegDomain} <: AbstractFEMM
    integdomain::ID # domain data
    mcsys::CSys # updater of the material orientation matrix
end

"""
    FEMMBase(integdomain::IntegDomain{S, F}) where {S<:AbstractFESet, F<:Function}

Construct with the default orientation matrix (identity).
"""
function FEMMBase(integdomain::ID) where {ID<:IntegDomain}
    return FEMMBase(integdomain, CSys(manifdim(integdomain.fes)))
end

"""
    finite_elements(self::AbstractFEMM)

Retrieve the finite element set for this FEMM to work on.
"""
function finite_elements(self::AbstractFEMM)
    self.integdomain.fes
end

"""
    associategeometry!(self::AbstractFEMM, geom::NodalField{FT}) where {FT}

Associate geometry field with the FEMM.

There may be operations that could benefit from pre-computations
that involve a geometry field. If so, associating the geometry
field gives the FEMM a chance to save on repeated computations.

Geometry field is normally passed into any routine that evaluates some
forms (integrals) over the mesh.  Whenever the geometry passed into a
routine is not consistent with the one for which `associategeometry!()`
was called before, `associategeometry!()` needs to be called with
the new geometry field.
"""
function associategeometry!(self::AbstractFEMM, geom::NodalField{FT}) where {FT}
    return self # default is no-op
end

"""
    inspectintegpoints(
        self::FEMM,
        geom::NodalField{FT},
        felist::AbstractVector{IT},
        inspector::F,
        idat,
        quantity = :Cauchy;
        context...,
    ) where {FEMM<:AbstractFEMM, FT, IT, F<:Function}

Inspect integration points.
"""
function inspectintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    felist::AbstractVector{IT},
    inspector::F,
    idat,
    quantity = :Cauchy;
    context...,
) where {FEMM<:AbstractFEMM, FT, IT, F<:Function}
    return idat # default is no-op
end

"""
    integratefieldfunction(
        self::AbstractFEMM,
        geom::NodalField{FT},
        afield::FL,
        fh::F,
        initial::R;
        m = -1,
    ) where {FT, FL<:NodalField{T}, T, F<:Function, R}

Integrate a nodal-field function over the discrete manifold.

`afield` = NODAL field to supply the values
`fh` = function taking position and the field value as arguments, returning
    value of type `T`.
`m` = dimension of the manifold over which to integrate; `m < 0` means that
    the dimension is controlled by the manifold dimension of the elements.

Returns value of type `R`, which is initialized by `initial`.
"""
function integratefieldfunction(
    self::AbstractFEMM,
    geom::NodalField{FT},
    afield::FL,
    fh::F,
    initial::R;
    m = -1,
) where {FT, T, FL<:NodalField{T}, F<:Function, R}
    fes = finite_elements(self)  # finite elements
    # Constants
    nfes = count(fes) # number of finite elements in the set
    ndn = ndofs(afield) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    a = fill(zero(T), nne, ndn) # array of field DOFS-- used as a buffer
    ecoords = fill(zero(FT), nne, ndofs(geom)) # array of field DOFS-- used as a buffer
    loc = fill(zero(FT), 1, sdim) # quadrature point location -- used as a buffer
    val = fill(zero(T), 1, ndn) # field value at the point -- used as a buffer
    J = fill(zero(FT), sdim, mdim) # Jacobian matrix -- used as a buffer
    if m >= 0
        # Either the manifold dimension was supplied
    else
        m = mdim# ...Or it is implied
    end
    result = initial           # initial value for the result
    for i in eachindex(fes) #Now loop over all fes in the block
        gathervalues_asmat!(afield, a, fes.conn[i])# retrieve element dofs
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        for j in 1:npts #Loop over all integration points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            my_At_mul_B!(val, Ns[j], a)# Field value at the quadrature point
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[i], Ns[j], m)
            result = result + fh(loc, val) * Jac * w[j]
        end
    end
    return result
end

"""
    integratefieldfunction(
        self::AbstractFEMM,
        geom::NodalField{FT},
        afield::FL,
        fh::F,
        initial::R;
        m = -1,
    ) where {FT, FL<:ElementalField{T}, T, F<:Function, R}

Integrate a elemental-field function over the discrete manifold.

`afield` = ELEMENTAL field to supply the values
`fh` = function taking position and the field value as arguments,
    returning value of type `T`.
`m` = dimension of the manifold over which to integrate; `m < 0` means that
    the dimension is controlled by the manifold dimension of the elements.

Returns value of type `T`, which is initialized by `initial`.

Returns value of type `R`, which is initialized by `initial`.
"""
function integratefieldfunction(
    self::AbstractFEMM,
    geom::NodalField{FT},
    afield::FL,
    fh::F,
    initial::R;
    m = -1,
) where {FT, T, FL<:ElementalField{T}, F<:Function, R}
    fes = finite_elements(self)  # finite elements
    # Constants
    nfes = count(fes) # number of finite elements in the set
    ndn = ndofs(afield) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    a = fill(zero(FT), nne, ndn) # array of field DOFS-- used as a buffer
    ecoords = fill(zero(FT), nne, ndofs(geom)) # array of field DOFS-- used as a buffer
    loc = fill(zero(FT), 1, sdim) # quadrature point location -- used as a buffer
    J = fill(zero(FT), sdim, mdim) # Jacobian matrix -- used as a buffer
    if m >= 0
        # Either the manifold dimension was supplied
    else
        m = mdim# ...Or it is implied
    end
    result = initial           # initial value for the result
    for i in eachindex(fes) #Now loop over all fes in the block
        gathervalues_asmat!(afield, a, [i])# retrieve element dofs
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        for j in 1:npts #Loop over all integration points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[i], Ns[j], m)
            result = result + fh(loc, a) * Jac * w[j]
        end
    end
    return result
end

"""
    integratefunction(
        self::AbstractFEMM,
        geom::NodalField{FT},
        fh::F,
        m = -1,
    ) where {FT, F<:Function}

Integrate a function over the discrete manifold.

Integrate some scalar function over the geometric cells. The function takes a
single argument, the position vector.

When the scalar function returns just +1 (such as `(x) ->  1.0`), the result
measures the volume (number of points, length, area, 3-D volume, according to
the manifold dimension). When the function returns the mass density, the
method measures the mass, when the function returns the x-coordinate equal
measure the static moment with respect to the y- axis, and so on.

# Example:
Compute the volume of the mesh and then its center of gravity:
```
V = integratefunction(femm, geom, (x) ->  1.0)
Sx = integratefunction(femm, geom, (x) ->  x[1])
Sy = integratefunction(femm, geom, (x) ->  x[2])
Sz = integratefunction(femm, geom, (x) ->  x[3])
CG = vec([Sx Sy Sz]/V)
```
Compute a moment of inertia of the mesh relative to the origin:
```
Ixx = integratefunction(femm, geom, (x) ->  x[2]^2 + x[3]^2)
```
"""
function integratefunction(
    self::AbstractFEMM,
    geom::NodalField{FT},
    fh::F,
    m = -1,
) where {FT, F<:Function}
    fes = finite_elements(self)
    if m < 0
        m = manifdim(fes)  # native  manifold dimension
    end
    # Constants
    nfes = count(fes) # number of finite elements in the set
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    ecoords = fill(zero(FT), nne, ndofs(geom)) # array of field DOFS-- used as a buffer
    loc = fill(zero(FT), 1, sdim) # quadrature point location -- used as a buffer
    J = fill(zero(FT), sdim, mdim) # Jacobian matrix -- used as a buffer
    result = 0.0# Initialize the result
    for i in eachindex(fes)  # Now loop over all fes in the set
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        for j in 1:npts #Loop over all integration points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[i], Ns[j], m)
            result = result + fh(vec(loc)) * Jac * w[j]
        end
    end
    return result
end

"""
    transferfield!(
        ff::F,
        fensf::FENodeSet{FT},
        fesf::AbstractFESet,
        fc::F,
        fensc::FENodeSet{FT},
        fesc::AbstractFESet,
        geometricaltolerance::FT;
        parametrictolerance::FT = 0.01,
    ) where {FT<:Number, F<:NodalField{T}, T}

Transfer a nodal field from a coarse mesh to a finer one.

# Arguments
- `ff` = the fine-mesh field (modified and also returned)
- `fensf` = finite element node set for the fine-mesh
- `fc` = the coarse-mesh field
- `fensc` = finite element node set for the fine-mesh,
- `fesc` = finite element set for the coarse mesh
- `geometricaltolerance` = tolerance in physical space for searches of the
  adjacent nodes
- `parametrictolerance` = tolerance in parametric space for for check whether
  node is inside an element

# Output
Nodal field `ff` transferred to the fine mesh is output.
"""
function transferfield!(
    ff::F,
    fensf::FENodeSet{FT},
    fesf::AbstractFESet,
    fc::F,
    fensc::FENodeSet{FT},
    fesc::AbstractFESet,
    geometricaltolerance::FT;
    parametrictolerance::FT = 0.01,
) where {FT<:Number, T, F<:NodalField{T}}
    fill!(ff.values, Inf) # the "infinity" value indicates a missed node
    @assert count(fensf) == nents(ff)
    parametrictol = 0.01
    nodebox = initbox!([], vec(fensc.xyz[1, :]))
    # Find out how many partitions of the nodes on the fine mesh we should use
    npartitions = max(2, Int(round(count(fensf) / 1000)))
    # Partition the nodes of the fine mesh
    npf = nodepartitioning(fensf, npartitions)
    partitionnumbers = unique(npf)
    npartitions = length(partitionnumbers)
    # Go through all the partitions
    for p in 1:npartitions
        pnl = findall(x -> x == partitionnumbers[p], npf) # subset of fine-mesh nodes
        # Find the bounding box
        subbox = boundingbox(fensf.xyz[pnl, :])
        tol = 2 * geometricaltolerance # increase the box a bit
        # Construct a sub mesh of the coarse mesh that covers the nodes from this partition
        sublist = selectelem(fensc, fesc, overlappingbox = subbox, inflate = tol)
        if !isempty(sublist) # there are some finite elements to work with
            fescsub = subset(fesc, sublist)
            connected = findunconnnodes(fensc, fescsub)
            fenscsub, newnumber = compactnodes(fensc, connected) # nodes of the sub mesh
            fescsub = renumberconn!(fescsub, newnumber) # elements of the sub mesh
            present = findall(x -> x > 0, newnumber)
            fcsub = NodalField(fc.values[present, :]) # reduce the coarse-mesh field to the sub mesh
            # Now we can find the values at the nodes of the subset of the fine mesh
            # working only with the sub mesh of the coarse mesh
            for i in pnl # for all nodes in the subset
                nl = vselect(fenscsub.xyz; nearestto = fensf.xyz[i, :])
                # For each node in the fine field try to find one in  the coarse field
                if !isempty(nl) &&
                   norm(fensf.xyz[i, :] - fenscsub.xyz[nl[1], :]) < geometricaltolerance
                    # These nodes correspond in the  refined  and coarse mesh
                    ff.values[i, :] = fcsub.values[nl[1], :]
                else
                    # Obviously, some nodes in the fine mesh are not located "at" the
                    # location of a coarse-mesh node. Then we have to search which
                    # element they fall into.
                    nodebox = initbox!(nodebox, vec(fensf.xyz[i, :]))
                    nodebox = inflatebox!(nodebox, geometricaltolerance)
                    el = selectelem(fenscsub, fescsub; overlappingbox = nodebox)
                    for e in el
                        c = [k for k in fescsub.conn[e]] #c = view(fescsub.conn, e, :)
                        pc, success = map2parametric(
                            fescsub,
                            fenscsub.xyz[c, :],
                            vec(fensf.xyz[i, :]);
                            tolerance = 0.000001,
                            maxiter = 7,
                        )
                        @assert success # this shouldn't be tripped; normally we succeed
                        if inparametric(fescsub, pc; tolerance = parametrictolerance) # coarse mesh element encloses the node
                            N = bfun(fescsub, pc)
                            ff.values[i, :] = transpose(N) * fcsub.values[c, :]
                            break
                        end
                    end
                end
            end # for i in pnl # for all nodes in the subset
        end # if !isempty(sublist)
    end # for p 

    # Check that we haven't missed any node connected to some finite elements.
    cnl = connectednodes(fesf)
    for i in cnl
        if any(v -> v == Inf, ff.values[i, :])
            println("fensf.xyz[$i, :] = $(fensf.xyz[i, :])")
            error("Missed node in transfer")
        end
    end

    return ff
end

"""
    transferfield!(
        ff::F,
        fensf::FENodeSet{FT},
        fesf::AbstractFESet,
        fc::F,
        fensc::FENodeSet{FT},
        fesc::AbstractFESet,
        geometricaltolerance::FT;
        parametrictolerance::FT = 0.01,
    ) where {FT<:Number, F<:ElementalField{T}, T}

Transfer an elemental field from a coarse mesh to a finer one.

# Arguments
- `ff` = the fine-mesh field (modified and also returned)
- `fensf` = finite element node set for the fine-mesh
- `fc` = the coarse-mesh field
- `fensc` = finite element node set for the fine-mesh,
- `fesc` = finite element set for the coarse mesh
- `tolerance` = tolerance in physical space for searches of the adjacent nodes

# Output
Elemental field `ff` transferred to the fine mesh is output.
"""
function transferfield!(
    ff::F,
    fensf::FENodeSet{FT},
    fesf::AbstractFESet,
    fc::F,
    fensc::FENodeSet{FT},
    fesc::AbstractFESet,
    geometricaltolerance::FT;
    parametrictolerance::FT = 0.01,
) where {FT<:Number, T, F<:ElementalField{T}}
    @assert count(fesf) == nents(ff)
    nodebox = initbox!([], vec(fensc.xyz[1, :]))
    centroidpc = centroidparametric(fesf)
    N = bfun(fesf, centroidpc)
    NT = transpose(N)
    for i in eachindex(fesf) # For all finite elements in the fine mesh
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        nodebox = initbox!(nodebox, vec(centroid))
        nodebox = inflatebox!(nodebox, geometricaltolerance)
        el = selectelem(fensc, fesc; overlappingbox = nodebox)
        foundone = false
        for e in el
            c = [k for k in fesc.conn[e]]
            pc, success = map2parametric(
                fesc,
                fensc.xyz[c, :],
                vec(centroid);
                tolerance = 0.000001,
                maxiter = 9,
            )
            # if !success
            # println("pc = $(pc)")
            # N1 = bfun(fesf, pc)
            # p = transpose(N1) * fensc.xyz[view(fesc.conn, e, :), :]
            # println("p = $(p)")
            # println("centroid = $(centroid)")
            # end
            # @assert success # this shouldn't be tripped; normally we succeed
            if success && inparametric(fesc, pc; tolerance = 0.001) # coarse mesh element encloses the centroid
                ff.values[i, :] = fc.values[e, :]
                foundone = true
                break
            end
        end
        @assert foundone
    end
    return ff
end

function  _buffers_basic(self, geom::NodalField{FT}, P::NodalField{T}) where {FT<:Number, T<:Number}
    fes = self.integdomain.fes
    ndn = ndofs(P); # number of degrees of freedom per node
    nne =  nodesperelem(fes); # number of nodes per element
    sdim =  ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes);     # manifold dimension of the element
    elmdim = ndn*nne;          # dimension of the element matrix
    ecoords = fill(zero(FT), nne, ndofs(geom)); # array of element coordinates
    IT = eltype(P.dofnums)
    dofnums = fill(zero(IT), elmdim); # degree of freedom array -- used as a buffer
    loc = fill(zero(FT), 1, sdim); # quadrature point location -- used as a buffer
    J = fill(zero(FT), sdim, mdim); # Jacobian matrix -- used as a buffer
    gradN = fill(zero(FT), nne, mdim); # intermediate result -- used as a buffer
    return nne, ndn, ecoords, dofnums, loc, J, gradN
end

function _buffers_el(self, geom::NodalField{FT}, P::NodalField{T}) where {FT<:Number, T<:Number}
    fes = self.integdomain.fes
    ndn = ndofs(P); # number of degrees of freedom per node
    nne =  nodesperelem(fes); # number of nodes per element
    sdim =  ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes);     # manifold dimension of the element
    elmdim = ndn*nne;          # dimension of the element matrix
    elmat = fill(zero(T), elmdim, elmdim);# element matrix -- used as a buffer
    elvec = fill(zero(T), elmdim); # buffer
    elvecfix = fill(zero(T), elmdim); # buffer
    return elmdim, elmat, elvec, elvecfix
end

"""
    linform_dot(
        self::FEMM,
        assembler::A,
        geom::NodalField{FT},
        P::NodalField{T},
        fi::ForceIntensity,
        m,
    ) where {FEMM<:AbstractFEMM, A<:AbstractSysvecAssembler, FT<:Number, T}

Evaluate the linear form "dot".

```math
\\int_{V}  \\vartheta \\cdot f \\; \\mathrm{d} V
```
Here ``\\vartheta`` is the test function, ``f`` is a given function (data).
Both can be vectors. ``f`` is represented with `ForceIntensity`.

# Arguments
- `fi`=force intensity object
- `m`= manifold dimension, 1= curve, 2= surface, 3= volume. For body loads `m`
is set to 3, for tractions on the surface it is set to 2, and so on.
"""
function linform_dot(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    P::NodalField{T},
    fi::ForceIntensity,
    m,
) where {FEMM<:AbstractFEMM, A<:AbstractSysvecAssembler, FT<:Number, T}
    fes = finite_elements(self)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    # Prepare some buffers:
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buffers_basic(self, geom, P)
    elmdim, elmat, elvec, elvecfix = _buffers_el(self, geom, P)
    startassembly!(assembler, P.nfreedofs)
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elvec, T(0.0))
        for j in 1:npts
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[i], Ns[j], m)
            force = updateforce!(fi, loc, J, fes.label[i]) # retrieve the applied load
            Factor = (Jac * w[j])
            NkxF = zero(T)
            rx = 1
            for kx in 1:nne # all the nodes
                NkxF = Ns[j][kx] * Factor
                for mx in 1:ndn   # all the degrees of freedom
                    elvec[rx] = elvec[rx] + NkxF * force[mx]
                    rx = rx + 1    # next component of the vector
                end
            end
        end
        gatherdofnums!(P, dofnums, fes.conn[i])
        assemble!(assembler, elvec, dofnums)
    end
    F = makevector!(assembler)
    return F
end


function linform_dot(
    self::FEMM,
    geom::NodalField{FT},
    P::NodalField{T},
    fi::ForceIntensity,
    m,
) where {FEMM<:AbstractFEMM, FT<:Number, T}
    assembler = SysvecAssembler(0.0 * P.values[1])#T(0.0))
    return distribloads(self, assembler, geom, P, fi, m)
end

# Alias
const distribloads  = linform_dot

"""
    connectionmatrix(self::FEMM, nnodes) where {FEMM<:AbstractFEMM}

Compute the connection matrix.

The matrix has a nonzero in all the rows and columns which correspond to nodes
connected by some finite element.
"""
function connectionmatrix(self::FEMM, nnodes) where {FEMM<:AbstractFEMM}
    fes = finite_elements(self)
    nfes = length(fes.conn)
    nconns = nodesperelem(fes)
    N = nfes * nconns * nconns
    IT = eltype(fes.conn[1])
    rb = IT[]
    sizehint!(rb, N)
    cb = IT[]
    sizehint!(cb, N)
    vb = ones(IT, N)
    for j in 1:nfes
        for k in 1:nconns
            append!(rb, fes.conn[j])
            @inbounds for m in 1:nconns
                push!(cb, fes.conn[j][k])
            end
        end
    end
    return sparse(rb, cb, vb, nnodes, nnodes)
end

"""
    dualconnectionmatrix(
        self::FEMM,
        fens::FENodeSet,
        minnodes = 1,
    ) where {FEMM<:AbstractFEMM}

Compute the dual connection matrix.

The matrix has a nonzero in all the rows and columns which correspond to
elements connected by some finite element nodes.

- `minnodes`: minimum number of nodes that the elements needs to share in order
  to be neighbors (default 1)
"""
function dualconnectionmatrix(
    self::FEMM,
    fens::FENodeSet,
    minnodes = 1,
) where {FEMM<:AbstractFEMM}
    fes = finite_elements(self)
    nfes = length(fes.conn)
    nconns = nodesperelem(fes)
    N = nfes * nconns * nconns
    IT = eltype(fes.conn[1])
    rb = IT[]
    sizehint!(rb, N)
    cb = IT[]
    sizehint!(cb, N)
    m = FENodeToFEMap(fes.conn, count(fens))
    for j in eachindex(m.map)
        for i in eachindex(m.map[j])
            append!(rb, m.map[j])
            @inbounds for k in eachindex(m.map[j])
                push!(cb, m.map[j][i])
            end
        end
    end
    vb = ones(IT, length(rb))
    C = sparse(rb, cb, vb, nfes, nfes)
    I, J, V = findnz(C)
    ix = findall(x -> x >= minnodes, V)
    return sparse(I[ix], J[ix], V[ix], nfes, nfes)
end


struct InverseDistanceInspectorData{IT, FT}
    component::Vector{IT}
    d::Vector{FT} # nodesperelem(integdomain.fes)
    sum_inv_dist::Vector{FT} # nnodes(geom)
    sum_quant_inv_dist::Array{FT,2} # nnodes(geom) x length(component)
end


# This is an inverse-distance interpolation inspector.
function _idi_inspector(idat, elnum, conn, xe, out, xq)
    # xe = coordinates of the nodes of the element
    # xq = coordinate of the quadrature point
    mindn = Inf
    for jjj in axes(xe, 1)
        idat.d[jjj] = sum((vec(xe[jjj, :]) - vec(xq)) .^ 2)
        if (idat.d[jjj] > 0.0)
            mindn = min(mindn, idat.d[jjj])
        end
    end
    mindn = mindn / 1.0e9
    for jjj in eachindex(idat.d)
        invdjjj = 1.0 / (idat.d[jjj] + mindn)
        quant = out[idat.component]
        for kkk in eachindex(quant)
            idat.sum_quant_inv_dist[conn[jjj], kkk] += invdjjj * quant[kkk]
        end
        idat.sum_inv_dist[conn[jjj]] += invdjjj
    end
    return idat
end


struct AveragingInspectorData{IT, FT}
    component::Vector{IT}
    d::Vector{FT} # nodesperelem(fes)
    ncontrib::Vector{IT} # nnodes(geom)
    sum_quant::Array{FT,2} # nnodes(geom) x length(component)
end

# This is a simple nodal averaging inspector.
# The quadrature point is assumed contribute only to the node that is nearest to
# it: The idea is that the inspected element  will produce outputs  at the
# locations of the nodes.
function _avg_inspector(idat, elnum, conn, xe, out, xq)
    # xe = coordinates of the nodes of the element
    # xq = coordinate of the quadrature point
    for jjj in axes(xe, 1)
        idat.d[jjj] = sum((vec(xe[jjj, :]) - vec(xq)) .^ 2)
    end
    # Find  the node nearest to the quadrature point
    (minval, ix) = findmin(idat.d)
    quant = out[idat.component]
    for kkk in eachindex(quant)
        idat.sum_quant[conn[ix], kkk] += quant[kkk]
    end
    idat.ncontrib[conn[ix]] += 1
    return idat
end

"""
    fieldfromintegpoints(
        self::FEMM,
        geom::NodalField{FT},
        u::NodalField{T},
        dT::NodalField{FT},
        quantity::Symbol,
        component::FIntVec;
        context...,
    ) where {FEMM<:AbstractFEMM, FT<:Number, T<:Number}

Construct nodal field from integration points.

# Arguments
- `geom`     - reference geometry field
- `u`        - displacement field
- `dT`       - temperature difference field
- `quantity`   - this is what you would assign to the 'quantity' argument
           of the material update!() method.
- `component`- component of the 'quantity' array: see the material update()
           method.
Keyword arguments
- `nodevalmethod` = `:invdistance` (the default) or `:averaging`;
- `reportat` = at which point should the  element quantities be reported?
    This argument is interpreted inside the `inspectintegpoints()` method.

# Output
- the new field that can be used to map values to colors and so on
"""
function fieldfromintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    dT::NodalField{FT},
    quantity::Symbol,
    component::AbstractVector{IT};
    context...,
) where {FEMM<:AbstractFEMM, FT<:Number, T<:Number, IT<:Integer}
    fes = finite_elements(self)
    # Constants
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    nodevalmethod = :invdistance
    reportat = :default
    for apair in pairs(context)
        sy, val = apair
        if sy == :nodevalmethod
            nodevalmethod = val
        end
        if sy == :reportat
            reportat = val
        end
    end
    if nodevalmethod == :averaging
        # Container of intermediate results
        idat = AveragingInspectorData(
            component,
            zeros(FT, nne),
            zeros(Int, nnodes(geom)),
            zeros(FT, nnodes(geom), length(component)),
        )
        # Loop over cells to interpolate to nodes
        idat = inspectintegpoints(
            self,
            geom,
            u,
            dT,
            collect(eachindex(fes)),
            _avg_inspector,
            idat,
            quantity;
            context...,
        )
        # The data for the field to be constructed is initialized
        nvals = zeros(FT, nnodes(geom), length(component))
        # compute the data array
        for kkk in axes(nvals, 1)
            for j in eachindex(component)
                if (idat.ncontrib[kkk] > 0)
                    nvals[kkk, j] = idat.sum_quant[kkk, j] / idat.ncontrib[kkk]
                end
            end
        end
    else # inverse distance
        @assert (reportat == :default) || (reportat == :meanonly) "Inverse-distance interpolation requires :meanonly"
        # Container of intermediate results
        idat = InverseDistanceInspectorData(
            component,
            zeros(FT, nne),
            zeros(FT, nnodes(geom)),
            zeros(FT, nnodes(geom), length(component)),
        )
        # Loop over cells to interpolate to nodes
        idat = inspectintegpoints(
            self,
            geom,
            u,
            dT,
            collect(eachindex(fes)),
            _idi_inspector,
            idat,
            quantity;
            context...,
        )
        # The data for the field to be constructed is initialized
        nvals = zeros(FT, nnodes(geom), length(component))
        # compute the data array
        for kkk in axes(nvals, 1)
            for j in eachindex(component)
                if (idat.sum_inv_dist[kkk] > 0.0)
                    nvals[kkk, j] = idat.sum_quant_inv_dist[kkk, j] / idat.sum_inv_dist[kkk]
                end
            end
        end
    end
    # Make the field
    return NodalField(nvals)
end

function fieldfromintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    dT::NodalField{FT},
    quantity::Symbol,
    component::IT;
    context...,
) where {FEMM<:AbstractFEMM, FT<:Number, T<:Number, IT<:Integer}
    return fieldfromintegpoints(self, geom, u, dT, quantity, [component]; context...)
end

function fieldfromintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    quantity::Symbol,
    component::AbstractVector{IT};
    context...,
) where {FEMM<:AbstractFEMM, FT<:Number, T<:Number, IT<:Integer}
    dT = NodalField(zeros(FT, nnodes(geom), 1)) # zero difference in temperature
    return fieldfromintegpoints(self, geom, u, dT, quantity, component; context...)
end

function fieldfromintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    quantity::Symbol,
    component::IT;
    context...,
) where {FEMM<:AbstractFEMM, FT<:Number, T<:Number, IT<:Integer}
    dT = NodalField(zeros(FT, nnodes(geom), 1)) # zero difference in temperature
    return fieldfromintegpoints(self, geom, u, dT, quantity, [component]; context...)
end


struct MeanValueInspectorData{IT, FT}
    n_quant::Vector{IT}
    sum_quant_value::Array{FT,2}
end

"""
    elemfieldfromintegpoints(
        self::FEMM,
        geom::NodalField{FT},
        u::NodalField{T},
        dT::NodalField{FT},
        quantity::Symbol,
        component::FIntVec;
        context...,
    ) where {FEMM<:AbstractFEMM, FT<:Number, T<:Number}

Construct elemental field from integration points.

# Arguments
`geom`     - reference geometry field
`u`        - displacement field
`dT`       - temperature difference field
`quantity`   - this is what you would assign to the 'quantity' argument
           of the material update!() method.
`component`- component of the 'quantity' array: see the material update()
           method.

# Output
 - the new field that can be used to map values to colors and so on
"""
function elemfieldfromintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    dT::NodalField{FT},
    quantity::Symbol,
    component::AbstractVector{IT};
    context...,
) where {FEMM<:AbstractFEMM, FT<:Number, T<:Number, IT<:Integer}
    fes = finite_elements(self)
    # Constants
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    # Container of intermediate results
    idat = MeanValueInspectorData(
        zeros(IT, count(fes)),
        zeros(FT, count(fes), length(component)),
    )
    # This is an mean-value interpolation inspector. The mean of the
    # quadrature-point quantities is reported per element.
    function mv_inspector(idat, elnum, conn, xe, out, xq)
        # xe = coordinates of the nodes of the element
        # xq = coordinate of the quadrature point
        idat.n_quant[elnum] += 1
        quant = out[component]
        for kkk in eachindex(quant)
            idat.sum_quant_value[elnum, kkk] += quant[kkk]
        end
        return idat
    end
    # Loop over cells to interpolate to nodes
    idat = inspectintegpoints(
        self,
        geom,
        u,
        dT,
        collect(eachindex(fes)),
        mv_inspector,
        idat,
        quantity;
        context...,
    )
    # The data for the field to be constructed is initialized
    evals = zeros(FT, count(fes), length(component))
    # compute the data array
    for j in axes(evals, 1)
        for kkk in axes(evals, 2)
            evals[j, kkk] = idat.sum_quant_value[j, kkk] / idat.n_quant[j]
        end
    end
    # Make the field
    return ElementalField(evals)
end

function elemfieldfromintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    dT::NodalField{FT},
    quantity::Symbol,
    component::IT;
    context...,
) where {FEMM<:AbstractFEMM, FT<:Number, T<:Number, IT<:Integer}
    return elemfieldfromintegpoints(self, geom, u, dT, quantity, [component]; context...)
end

function elemfieldfromintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    quantity::Symbol,
    component::IT;
    context...,
) where {FEMM<:AbstractFEMM, FT<:Number, T<:Number, IT<:Integer}
    dT = NodalField(zeros(FT, nnodes(geom), 1)) # zero difference in temperature
    return elemfieldfromintegpoints(self, geom, u, dT, quantity, [component]; context...)
end

function elemfieldfromintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    quantity::Symbol,
    component::AbstractVector{IT};
    context...,
) where {FEMM<:AbstractFEMM, FT<:Number, T<:Number, IT<:Integer}
    dT = NodalField(zeros(FT, nnodes(geom), 1)) # zero difference in temperature
    return elemfieldfromintegpoints(self, geom, u, dT, quantity, component; context...)
end


function _buffers(
    self::FEMM,
    geom::NodalField{FT},
    afield::NodalField{T},
) where {FEMM<:AbstractFEMM, FT, T}
    # Constants
    fes = finite_elements(self)
    nfes = count(fes) # number of finite elements in the set
    ndn = ndofs(afield) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)   # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    Kedim = ndn * nne      # dimension of the element matrix
    elmat = fill(zero(FT), Kedim, Kedim) # buffer
    ecoords = fill(zero(FT), nne, ndofs(geom)) # array of Element coordinates
    dofnums = fill(zero(eltype(afield.dofnums)), Kedim) # buffer
    loc = fill(zero(FT), 1, sdim) # buffer
    J = fill(zero(FT), sdim, mdim) # buffer
    gradN = fill(zero(FT), nne, mdim) # buffer
    return ecoords, dofnums, loc, J, gradN, elmat
end

"""
    innerproduct(
        self::FEMM,
        assembler::A,
        geom::NodalField{FT},
        afield::NodalField{T},
    ) where {FEMM<:AbstractFEMM, A<:AbstractSysmatAssembler, FT, T}

Compute the inner-product (Gram) matrix.
"""
function innerproduct(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    afield::NodalField{T},
) where {FEMM<:AbstractFEMM, A<:AbstractSysmatAssembler, FT, T}
    fes = finite_elements(self)
    ecoords, dofnums, loc, J, gradN, elmat = _buffers(self, geom, afield)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    NexpTNexp = Matrix{eltype(Ns[1])}[]# basis f. matrix -- buffer
    ndn = ndofs(afield)
    Indn = [i == j ? one(FT) : zero(FT) for i in 1:ndn, j in 1:ndn] # "identity"
    for j in 1:npts # This quantity is the same for all quadrature points
        Nexp = fill(zero(FT), ndn, size(elmat, 1))
        for l1 in 1:nodesperelem(fes)
            Nexp[1:ndn, (l1-1)*ndn+1:(l1)*ndn] = Indn * Ns[j][l1]
        end
        push!(NexpTNexp, Nexp' * Nexp)
    end
    startassembly!(
        assembler,
        size(elmat, 1),
        size(elmat, 2),
        count(fes),
        afield.nfreedofs,
        afield.nfreedofs,
    )
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elmat, 0.0) # Initialize element matrix
        for j in 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            thefactor::FT = (Jac * w[j])
            elmat .+= NexpTNexp[j] * thefactor
        end # Loop over quadrature points
        gatherdofnums!(afield, dofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function innerproduct(
    self::FEMM,
    geom::NodalField{FT},
    afield::NodalField{T},
) where {FEMM<:AbstractFEMM, FT, T}
    assembler = SysmatAssemblerSparseSymm()
    return innerproduct(self, assembler, geom, afield)
end


"""
    field_elem_to_nodal!(
        self::AbstractFEMM,
        geom::NodalField{FT},
        ef::EFL,
        nf::NFL;
        kind = :weighted_average,
    ) where {FT, T<:Number, EFL<:ElementalField{T}, NFL<:NodalField{T}}

Make a nodal field  from an elemental field over the discrete manifold.

`ef` = ELEMENTAL field to supply the values
`nf` = NODAL field to receive the values
`kind` = default is `:weighted_average`; other options: `:max`

Returns `nf`.
"""
function field_elem_to_nodal!(
    self::AbstractFEMM,
    geom::NodalField{FT},
    ef::EFL,
    nf::NFL;
    kind = :weighted_average,
) where {FT, T<:Number, EFL<:ElementalField{T}, NFL<:NodalField{T}}
    if kind == :max
        return _field_elem_to_nodal_max!(self, geom, ef, nf)
    else #:weighted_average
        return _field_elem_to_nodal_weighted_average!(self, geom, ef, nf)
    end
end

function _field_elem_to_nodal_weighted_average!(
    self::AbstractFEMM,
    geom::NodalField{FT},
    ef::EFL,
    nf::NFL,
) where {FT, T<:Number, EFL<:ElementalField{T}, NFL<:NodalField{T}}
    fes = finite_elements(self)  # finite elements
    # Dimensions
    nfes = count(fes) # number of finite elements in the set
    ndn = ndofs(nf) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    a = fill(zero(T), nne, ndn) # array of field DOFS-- used as a buffer
    ecoords = fill(zero(FT), nne, ndofs(geom)) # array of field DOFS-- used as a buffer
    loc = fill(zero(FT), 1, sdim) # quadrature point location -- used as a buffer
    J = fill(zero(FT), sdim, mdim) # Jacobian matrix -- used as a buffer
    nvolums = fill(zero(FT), nents(nf))
    # initial value for the result
    nf.values .= zero(T)
    for i in eachindex(fes) #Now loop over all fes in the block
        ev = ef.values[i, :]
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        for j in 1:npts #Loop over all integration points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[i], Ns[j], mdim)
            for k in eachindex(fes.conn[i])
                g = fes.conn[i][k]
                nf.values[g, :] .+= ev * Jac * w[j]
                nvolums[g] += Jac * w[j]
            end
        end
    end
    for g in 1:nents(nf)
        nf.values[g, :] ./= nvolums[g]
    end
    return nf
end

function _field_elem_to_nodal_max!(
    self::AbstractFEMM,
    geom::NodalField{FT},
    ef::EFL,
    nf::NFL,
) where {FT, T<:Number, EFL<:ElementalField{T}, NFL<:NodalField{T}}
    fes = finite_elements(self)  # finite elements
    nf.values .= zero(T) - Inf
    for i in eachindex(fes) #Now loop over all fes in the block
        ev = ef.values[i, :]
        for k in eachindex(fes.conn[i])
            g = fes.conn[i][k]
            nf.values[g, :] .= max.(ev, nf.values[g, :])
        end
    end
    return nf
end

"""
    field_nodal_to_elem!(
        self::AbstractFEMM,
        geom::NodalField{FT},
        nf::NFL,
        ef::EFL;
        kind = :weighted_average,
    ) where {FT<:Number, T, EFL<:ElementalField{T}, NFL<:NodalField{T}}

Make an elemental field  from a nodal field over the discrete manifold.

`nf` = NODAL field to supply the values
`ef` = ELEMENTAL field to receive the values
`kind` = default is `:weighted_average`; other options: `:max`

Returns `ef`.
"""
function field_nodal_to_elem!(
    self::AbstractFEMM,
    geom::NodalField{FT},
    nf::NFL,
    ef::EFL;
    kind = :weighted_average,
) where {FT<:Number, T, EFL<:ElementalField{T}, NFL<:NodalField{T}}
    if kind == :max
        return _field_nodal_to_elem_max!(self, geom, nf, ef)
    else #:weighted_average
        return _field_nodal_to_elem_weighted_average!(self, geom, nf, ef)
    end
end

function _field_nodal_to_elem_weighted_average!(
    self::AbstractFEMM,
    geom::NodalField{FT},
    nf::NFL,
    ef::EFL,
) where {FT<:Number, T, EFL<:ElementalField{T}, NFL<:NodalField{T}}
    fes = finite_elements(self)  # finite elements
    # Dimensions
    nfes = count(fes) # number of finite elements in the set
    ndn = ndofs(ef) # number of degrees of freedom per element
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    a = fill(zero(T), nne, ndn) # array of field DOFS-- used as a buffer
    ev = fill(zero(T), ndn) # array of field DOFS-- used as a buffer
    ecoords = fill(zero(FT), nne, ndofs(geom)) # array of field DOFS-- used as a buffer
    loc = fill(zero(FT), 1, sdim) # quadrature point location -- used as a buffer
    J = fill(zero(FT), sdim, mdim) # Jacobian matrix -- used as a buffer
    nvolums = fill(zero(FT), nents(nf))
    # initial value for the result
    ef.values .= zero(T)
    for i in eachindex(fes) #Now loop over all fes in the block
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        for j in 1:npts #Loop over all integration points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[i], Ns[j], mdim)
            for k in eachindex(fes.conn[i])
                g = fes.conn[i][k]
                nvolums[g] += Ns[j][k] * Jac * w[j]
            end
        end
    end
    for i in eachindex(fes) #Now loop over all fes in the block
        gathervalues_asmat!(nf, a, fes.conn[i])# retrieve element dofs
        ev .= zero(T)
        evol = 0.0
        for k in eachindex(fes.conn[i])
            g = fes.conn[i][k]
            ev .+= a[k] * nvolums[g]
            evol += nvolums[g]
        end
        ef.values[i, :] .= ev / evol
    end
    return ef
end

function _field_nodal_to_elem_max!(
    self::AbstractFEMM,
    geom::NodalField{FT},
    nf::NFL,
    ef::EFL,
) where {FT<:Number, T, EFL<:ElementalField{T}, NFL<:NodalField{T}}
    fes = finite_elements(self)  # finite elements
    # Dimensions
    nfes = count(fes) # number of finite elements in the set
    ndn = ndofs(ef) # number of degrees of freedom per element
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    a = fill(zero(T), nne, ndn) # array of field DOFS-- used as a buffer
    ev = fill(zero(T), ndn) # array of field DOFS-- used as a buffer
    # initial value for the result
    ef.values .= zero(T) - Inf
    for i in eachindex(fes) #Now loop over all fes in the block
        gathervalues_asmat!(nf, a, fes.conn[i])# retrieve element dofs
        ev .= zero(T) - Inf
        for k in eachindex(fes.conn[i])
            g = fes.conn[i][k]
            ev .= max.(ev, a[k])
        end
        ef.values[i, :] = ev
    end
    return ef
end


"""
    acousticmass(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}

Compute the acoustic mass matrix.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler
- `geom` = geometry field
- `P` = acoustic (perturbation) pressure field

Return a matrix.
"""
# function bilinear_diffusion(
#     self::FEMM,
#     assembler::A,
#     geom::NodalField{FT},
#     P::NodalField{TT}
#     ) where {FEMM<:AbstractFEMM, A<:AbstractSysvecAssembler, FT<:Number, TT<:Number}
#     fes = self.integdomain.fes
#     ecoords, dofnums, loc, J, gradN, elmat, elvec, elvecfix =   buffers(self, geom, P)
#     # Precompute basis f. values + basis f. gradients wrt parametric coor
#     npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
#     Jac = 0.0;
#     afactor = T(0.0);
#     startassembly!(assembler, size(elmat,1), size(elmat,2), count(fes),
#         P.nfreedofs, P.nfreedofs);
#     for i = 1:count(fes) # Loop over elements
#         gathervalues_asmat!(geom, ecoords, fes.conn[i]);
#         fill!(elmat, T(0.0)); # Initialize element matrix
#         for j = 1:npts # Loop over quadrature points
#             locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
#             Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
#             # gradient WRT global Cartesian coordinates
#             gradN!(fes, gradN, gradNparams[j], J);
#             afactor = (Jac*w[j]);
#             add_mggt_ut_only!(elmat, gradN, afactor)
#         end # Loop over quadrature points
#         complete_lt!(elmat)
#         gatherdofnums!(P, dofnums, fes.conn[i]);# retrieve degrees of freedom
#         assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
#     end # Loop over elements
#     return makematrix!(assembler);
# end

# function bilinear_diffusion(self::FEMMAcoust, geom::NodalField, P::NodalField{T}) where {T<:Number}
#     # Make the default assembler object.
#     assembler  =  SysmatAssemblerSparseSymm();
#     return acousticmass(self, assembler, geom, P);
# end


end # module



