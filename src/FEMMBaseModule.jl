"""
    FEMMBaseModule

Module for comments/base operations on interiors and boundaries of domains.
"""
module FEMMBaseModule

__precompile__(true)

using LinearAlgebra
using LinearAlgebra: mul!, Transpose
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
using SparseArrays: sparse, findnz
using LinearAlgebra: norm
using ..FENodeSetModule: FENodeSet
using ..FESetModule:
    AbstractFESet,
    manifdim,
    nodesperelem,
    subset,
    map2parametric,
    inparametric,
    centroidparametric,
    bfun,
    gradN!
using ..IntegDomainModule: IntegDomain, integrationdata, Jacobianmdim, Jacobianvolume
using ..CSysModule: CSys, updatecsmat!, csmat
using ..FieldModule: ndofs, nents, gatherdofnums!, gathervalues_asmat!, nalldofs, nfreedofs
using ..NodalFieldModule: NodalField, nnodes
using ..ElementalFieldModule: ElementalField, nelems
using ..ForceIntensityModule: ForceIntensity, updateforce!
using ..DataCacheModule: DataCache
using ..MatrixUtilityModule: locjac!, add_nnt_ut_only!
using ..MatrixUtilityModule:
    add_gkgt_ut_only!,
    complete_lt!,
    locjac!,
    add_nnt_ut_only!,
    mulCAtB!,
    mulCAB!,
    add_mggt_ut_only!,
    add_btdb_ut_only!
using ..BoxModule: initbox!, boundingbox, inflatebox!
using ..MeshModificationModule: nodepartitioning, compactnodes, renumberconn!
using ..MeshSelectionModule: selectelem, vselect, findunconnnodes, connectednodes
using ..AssemblyModule:
    AbstractSysvecAssembler,
    AbstractSysmatAssembler,
    SysmatAssemblerSparse,
    SysmatAssemblerSparseSymm,
    startassembly!,
    assemble!,
    makematrix!,
    makevector!,
    SysvecAssembler
using ..FENodeToFEMapModule: FENodeToFEMap
using ..DeforModelRedModule: AbstractDeforModelRed, blmat!, nstressstrain


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
) where {FEMM<:AbstractFEMM,FT,IT,F<:Function}
    return idat # default is no-op
end

"""
    integratefieldfunction(
        self::AbstractFEMM,
        geom::NodalField{FT},
        afield::FL,
        fh::F;
        initial::R,
        m = -1,
    ) where {FT, FL<:NodalField{T}, T, F<:Function, R}

Integrate a nodal-field function over the discrete manifold.

- `afield` = NODAL field to supply the values
- `fh` = function taking position and an array of field values for the element
  as arguments, returning value of type `T`. The rectangular array of field
  values has as many rows as there are nodes per element, and as many columns
  as there are degrees of freedom per node.
- `m` = dimension of the manifold over which to integrate; `m < 0` means that
  the dimension is controlled by the manifold dimension of the elements.

Returns value of type `R`, which is initialized by `initial`.
"""
function integratefieldfunction(
    self::AbstractFEMM,
    geom::NodalField{FT},
    afield::FL,
    fh::F;
    initial::R,
    m = -1,
) where {FT,T,FL<:NodalField{T},F<:Function,R}
    fes = finite_elements(self)
    if m < 0
        m = manifdim(fes)  # native  manifold dimension
    end
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buff_b(self, geom, afield)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    a = fill(zero(eltype(afield.values)), nne, ndn) # used as a buffer
    val = fill(zero(eltype(afield.values)), 1, ndn) # used as a buffer
    function fillcache!(
        cacheout::R,
        XYZ::VecOrMat{T},
        tangents::Matrix{T},
        feid::IT, qpid::IT
    ) where {R,T,IT}
        gathervalues_asmat!(afield, a, fes.conn[feid])# retrieve element dofs
        At_mul_B!(val, Ns[qpid], a)# Field value at the quadrature point
        cacheout = fh(XYZ, val)
        return cacheout
    end
    return ev_integrate(self, geom, DataCache(initial, fillcache!), initial, m)
end

"""
    integratefieldfunction(
        self::AbstractFEMM,
        geom::NodalField{FT},
        afield::FL,
        fh::F;
        initial::R = zero(FT),
        m = -1,
    ) where {FT, T, FL<:ElementalField{T}, F<:Function, R}

Integrate a elemental-field function over the discrete manifold.

- `afield` = ELEMENTAL field to supply the values
- `fh` = function taking position and an array of field values for the element
  as arguments, returning value of type `T`. The array of field values has one
  row and as many columns as there are degrees of freedom per element.
- `m` = dimension of the manifold over which to integrate; `m < 0` means that
  the dimension is controlled by the manifold dimension of the elements.

Returns value of type `R`, which is initialized by `initial`.
"""
function integratefieldfunction(
    self::AbstractFEMM,
    geom::NodalField{FT},
    afield::FL,
    fh::F;
    initial::R = zero(FT),
    m = -1,
) where {FT,T,FL<:ElementalField{T},F<:Function,R}
    fes = finite_elements(self)
    if m < 0
        m = manifdim(fes)  # native  manifold dimension
    end
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buff_b(self, geom, afield)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    a = fill(zero(eltype(afield.values)), nne, ndn) # used as a buffer
    function fillcache!(
        cacheout::R,
        XYZ::VecOrMat{T},
        tangents::Matrix{T},
        feid::IT, qpid::IT
    ) where {R,T,IT}
        gathervalues_asmat!(afield, a, (feid,))# retrieve element dofs
        cacheout = fh(XYZ, a)
        return cacheout
    end
    return ev_integrate(self, geom, DataCache(initial, fillcache!), initial, m)
end

"""
    integratefunction(
        self::AbstractFEMM,
        geom::NodalField{FT},
        fh::F;
        initial::R = zero(FT),
        m = -1,
    ) where {FT, F<:Function, R<:Number}

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
    fh::F;
    initial::R = zero(FT),
    m = -1,
) where {FT, F<:Function, R<:Number}
    fes = finite_elements(self)
    if m < 0
        m = manifdim(fes)  # native  manifold dimension
    end
    function fillcache!(
        cacheout::R,
        XYZ::VecOrMat{T},
        tangents::Matrix{T},
        feid::IT, qpid::IT
    ) where {R,T,IT}
        cacheout = fh(XYZ)
        return cacheout
    end
    return ev_integrate(self, geom, DataCache(initial, fillcache!), initial, m)
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
) where {FT<:Number,T,F<:NodalField{T}}
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
    for p = 1:npartitions
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
) where {FT<:Number,T,F<:ElementalField{T}}
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
    for j = 1:nfes
        for k = 1:nconns
            append!(rb, fes.conn[j])
            @inbounds for m = 1:nconns
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


struct InverseDistanceInspectorData{IT,FT}
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


struct AveragingInspectorData{IT,FT}
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
) where {FEMM<:AbstractFEMM,FT<:Number,T<:Number,IT<:Integer}
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
) where {FEMM<:AbstractFEMM,FT<:Number,T<:Number,IT<:Integer}
    return fieldfromintegpoints(self, geom, u, dT, quantity, [component]; context...)
end

function fieldfromintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    quantity::Symbol,
    component::AbstractVector{IT};
    context...,
) where {FEMM<:AbstractFEMM,FT<:Number,T<:Number,IT<:Integer}
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
) where {FEMM<:AbstractFEMM,FT<:Number,T<:Number,IT<:Integer}
    dT = NodalField(zeros(FT, nnodes(geom), 1)) # zero difference in temperature
    return fieldfromintegpoints(self, geom, u, dT, quantity, [component]; context...)
end


struct MeanValueInspectorData{IT,FT}
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
) where {FEMM<:AbstractFEMM,FT<:Number,T<:Number,IT<:Integer}
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
) where {FEMM<:AbstractFEMM,FT<:Number,T<:Number,IT<:Integer}
    return elemfieldfromintegpoints(self, geom, u, dT, quantity, [component]; context...)
end

function elemfieldfromintegpoints(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    quantity::Symbol,
    component::IT;
    context...,
) where {FEMM<:AbstractFEMM,FT<:Number,T<:Number,IT<:Integer}
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
) where {FEMM<:AbstractFEMM,FT<:Number,T<:Number,IT<:Integer}
    dT = NodalField(zeros(FT, nnodes(geom), 1)) # zero difference in temperature
    return elemfieldfromintegpoints(self, geom, u, dT, quantity, component; context...)
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
) where {FT,T<:Number,EFL<:ElementalField{T},NFL<:NodalField{T}}
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
) where {FT,T<:Number,EFL<:ElementalField{T},NFL<:NodalField{T}}
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
        for j = 1:npts #Loop over all integration points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[i], Ns[j], mdim)
            for k in eachindex(fes.conn[i])
                g = fes.conn[i][k]
                nf.values[g, :] .+= ev * Jac * w[j]
                nvolums[g] += Jac * w[j]
            end
        end
    end
    for g = 1:nents(nf)
        nf.values[g, :] ./= nvolums[g]
    end
    return nf
end

function _field_elem_to_nodal_max!(
    self::AbstractFEMM,
    geom::NodalField{FT},
    ef::EFL,
    nf::NFL,
) where {FT,T<:Number,EFL<:ElementalField{T},NFL<:NodalField{T}}
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
) where {FT<:Number,T,EFL<:ElementalField{T},NFL<:NodalField{T}}
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
) where {FT<:Number,T,EFL<:ElementalField{T},NFL<:NodalField{T}}
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
        for j = 1:npts #Loop over all integration points
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
) where {FT<:Number,T,EFL<:ElementalField{T},NFL<:NodalField{T}}
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


function _buff_b(self, geom, P)
    fes = finite_elements(self)
    ndn = ndofs(P) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    elmdim = ndn * nne          # dimension of the element matrix
    FT = eltype(geom.values)
    ecoords = fill(zero(FT), nne, ndofs(geom)) # array of element coordinates
    IT = eltype(P.dofnums)
    dofnums = fill(zero(IT), elmdim) # degree of freedom array -- used as a buffer
    loc = fill(zero(FT), 1, sdim) # quadrature point location -- used as a buffer
    J = fill(zero(FT), sdim, mdim) # Jacobian matrix -- used as a buffer
    gradN = fill(zero(FT), nne, mdim) # intermediate result -- used as a buffer
    return nne, ndn, ecoords, dofnums, loc, J, gradN
end

function _buff_e(self, geom, P, assembler)
    fes = finite_elements(self)
    ndn = ndofs(P) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    elmdim = ndn * nne          # dimension of the element matrix
    T = eltype(assembler) # type of the elementwise matrices, and of the global one as well
    elmat = fill(zero(T), elmdim, elmdim)# element matrix -- used as a buffer
    elvec = fill(zero(T), elmdim) # buffer
    return elmdim, elmat, elvec
end

"""
    ev_integrate(
            self::FEMM,
            geom::NodalField{FT},
            f::DC,
            initial::R,
            m,
    ) where {FEMM<:AbstractFEMM, FT<:Number, DC<:DataCache, R}

Compute the integral of a given function over a mesh domain.

```math
\\int_{\\Omega}  {f} \\; \\mathrm{d} \\Omega
```

Here ``{f}`` is a given function (data).  The data ``{f}`` is
represented with [`DataCache`](@ref).

# Arguments
- `self` = finite element machine;
- `geom` = geometry field;
- `f`= data cache, which is called to evaluate the integrand based on the
  location, the Jacobian matrix, the finite element identifier, and the
  quadrature point;
- `initial` = initial value of the integral,
- `m`= manifold dimension, 0= vertex (point), 1= curve, 2= surface, 3= volume.
  For body loads `m` is set to 3, for tractions on the surface it is set to 2,
  and so on.
"""
function ev_integrate(
    self::FEMM,
    geom::NodalField{FT},
    f::DC,
    initial::R,
    m,
) where {FEMM<:AbstractFEMM,FT<:Number,DC<:DataCache,R}
    fes = finite_elements(self)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buff_b(self, geom, geom)
    result = initial # Initialize the result
    for i in eachindex(fes)  # Now loop over all fes in the set
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        for j  in 1:npts #Loop over all integration points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[i], Ns[j], m)
            val = f(loc, J, i, j)
            result = result + val * Jac * w[j]
        end
    end
    return result
end

"""
    linform_dot(
        self::FEMM,
        assembler::A,
        geom::NodalField{FT},
        P::NodalField{T},
        f::DC,
        m,
    ) where {FEMM<:AbstractFEMM, A<:AbstractSysvecAssembler, FT<:Number, T, DC<:DataCache}

Compute the discrete vector implied by the linear form "dot".

```math
\\int_{V}  \\mathbf{w} \\cdot \\mathbf{f} \\; \\mathrm{d} V
```

Here ``\\mathbf{w}`` is the test function, ``\\mathbf{f}`` is a given function
(data). Both are assumed to be vectors,  even if they are of length 1,
representing scalars. The data ``\\mathbf{f}`` is represented with [`DataCache`]
(@ref).

# Arguments
- `self` = finite element machine;
- `assembler` = assembler of the global vector;
- `geom` = geometry field;
- `P` = nodal field to define the degree of freedom numbers;
- `f`= data cache, which is called to evaluate the integrand based on the
  location, the Jacobian matrix, the finite element identifier, and the
  quadrature point;
- `m`= manifold dimension, 0= vertex (point), 1= curve, 2= surface, 3= volume.
  For body loads `m` is set to 3, for tractions on the surface it is set to 2,
  and so on.
"""
function linform_dot(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    P::NodalField{T},
    f::DC,
    m,
) where {FEMM<:AbstractFEMM,A<:AbstractSysvecAssembler,FT<:Number,T,DC<:DataCache}
    fes = finite_elements(self)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    # Prepare some buffers:
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buff_b(self, geom, P)
    elmdim, elmat, elvec = _buff_e(self, geom, P, assembler)
    startassembly!(assembler, nalldofs(P))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elvec, T(0.0))
        for j  in  1:npts
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[i], Ns[j], m)
            force = f(loc, J, i, j) # retrieve the applied load
            Factor = (Jac * w[j])
            NkxF = zero(T)
            rx = 1
            for kx = 1:nne # all the nodes
                NkxF = Ns[j][kx] * Factor
                for mx = 1:ndn   # all the degrees of freedom
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
    f::DC,
    m,
) where {FEMM<:AbstractFEMM,FT<:Number,T,DC<:DataCache}
    assembler = SysvecAssembler(zero(eltype(P.values)))
    return linform_dot(self, assembler, geom, P, f, m)
end

"""
    distribloads(
        self::FEMM,
        assembler::A,
        geom::NodalField{FT},
        P::NodalField{T},
        fi::ForceIntensity,
        m,
    ) where {FEMM<:AbstractFEMM, A<:AbstractSysvecAssembler, FT<:Number, T}

Compute distributed loads vector.

# Arguments
- `fi`=force intensity object
- `m`= manifold dimension, 0= vertex (point), 1= curve, 2= surface, 3= volume.
  For body loads `m` is set to 3, for tractions on the surface it is set to 2,
  and so on.

The actual work is done by `linform_dot()`.
"""
function distribloads(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    P::NodalField{T},
    fi::ForceIntensity,
    m,
) where {FEMM<:AbstractFEMM,A<:AbstractSysvecAssembler,FT<:Number,T}
    return linform_dot(self, assembler, geom, P, fi._cache, m)
end

function distribloads(
    self::FEMM,
    geom::NodalField{FT},
    P::NodalField{T},
    fi::ForceIntensity,
    m,
) where {FEMM<:AbstractFEMM,FT<:Number,T}
    assembler = SysvecAssembler(zero(eltype(P.values)))
    return distribloads(self, assembler, geom, P, fi, m)
end

"""
    bilform_dot(
        self::FEMM,
        assembler::A,
        geom::NodalField{FT},
        u::NodalField{T},
        cf::DC
    ) where {FEMM<:AbstractFEMM, A<:AbstractSysmatAssembler, FT, T, DC<:DataCache}

Compute the sparse matrix implied by the bilinear form of the "dot" type.

```math
\\int_{\\Omega}  \\mathbf{w} \\cdot \\mathbf{c} \\cdot \\mathbf{u} \\; \\mathrm{d} \\Omega
```

Here ``\\mathbf{w}`` is the test function, ``\\mathbf{u}`` is the trial
function, ``\\mathbf{c}`` is a square matrix of coefficients; ``\\mathbf
{c}`` is computed by `cf`, which is a given function (data). Both trial and
test functions are assumed to be vectors(even if of length 1). `cf` is
represented with [`DataCache`](@ref), and needs to return a square matrix, with
dimension equal to the number of degrees of freedom per node in the `u` field.

The integral domain ``\\Omega`` can be the volume of the domain ``V`` (i.e. a
three dimensional integral), or a surface ``S`` (i.e. a two dimensional
integral), or a line domain ``L`` (i.e. a one dimensional integral).

# Arguments
- `self` = finite element machine;
- `assembler` = assembler of the global object;
- `geom` = geometry field;
- `u` = nodal field to define the degree of freedom numbers;
- `cf`= data cache, which is called to evaluate the coefficient ``c``, given the
  location of the integration point, the Jacobian matrix, and the finite
  element label.
- `m` = manifold dimension (default is 3).
"""
function bilform_dot(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    u::NodalField{T},
    cf::DC;
    m = 3,
) where {FEMM<:AbstractFEMM,A<:AbstractSysmatAssembler,FT,T,DC<:DataCache}
    fes = finite_elements(self)
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buff_b(self, geom, u)
    elmdim, elmat, elvec = _buff_e(self, geom, u, assembler)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    startassembly!(assembler, prod(size(elmat)) * count(fes), nalldofs(u), nalldofs(u))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elmat, 0.0) # Initialize element matrix
        for j  in  1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianmdim(self.integdomain, J, loc, fes.conn[i], Ns[j], m)
            c = cf(loc, J, i, j)
            for k = 1:nne, m = 1:nne
                factor = (Ns[j][k] * Ns[j][m] * Jac * w[j])
                for p  in  1:ndn, q  in  1:ndn
                    elmat[(k-1)*ndn+p, (m-1)*ndn+q] += factor * c[p, q]
                end
            end
        end # Loop over quadrature points
        gatherdofnums!(u, dofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function bilform_dot(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    cf::DC,
) where {FEMM<:AbstractFEMM,FT,T,DC<:DataCache}
    assembler = SysmatAssemblerSparseSymm()
    return bilform_dot(self, assembler, geom, u, cf)
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
) where {FEMM<:AbstractFEMM,A<:AbstractSysmatAssembler,FT,T}
    return bilform_dot(
        self,
        assembler,
        geom,
        afield,
        DataCache(LinearAlgebra.I(ndofs(afield))),
    )
end

function innerproduct(
    self::FEMM,
    geom::NodalField{FT},
    afield::NodalField{T},
) where {FEMM<:AbstractFEMM,FT,T}
    assembler = SysmatAssemblerSparseSymm()
    return innerproduct(self, assembler, geom, afield)
end

function _buff_d(self, geom, u)
    fes = finite_elements(self)
    nne = nodesperelem(fes) # number of nodes for element
    mdim = manifdim(fes) # manifold dimension of the element
    FT = eltype(geom.values)
    RmTJ = fill(zero(FT), mdim, mdim) # buffer
    T = eltype(u.values)
    c_gradNT = fill(zero(T), mdim, nne) # buffer
    return RmTJ, c_gradNT
end

"""
    bilform_diffusion(
        self::FEMM,
        assembler::A,
        geom::NodalField{FT},
        u::NodalField{T},
        cf::DC
    ) where {FEMM<:AbstractFEMM, A<:AbstractSysmatAssembler, FT, T, DC<:DataCache}

Compute the sparse matrix implied by the bilinear form of the "diffusion" type.

```math
\\int_{V}  \\nabla w \\cdot c \\cdot \\nabla u \\; \\mathrm{d} V
```

Here ``\\nabla w`` is the gradient of the scalar test function, ``\\nabla u`` is
the gradient of the scalar trial function, ``c`` is a square symmetric matrix
of coefficients or a scalar; ``c`` is computed by `cf`, which is a given
function (data). Both test and trial functions are assumed to be from the same
approximation space. `cf` is represented with [`DataCache`](@ref), and needs to
return a symmetric square matrix (to represent general anisotropic diffusion)
or a scalar (to represent isotropic diffusion).

The coefficient matrix ``c`` can be given in the so-called local material
coordinates: coordinates that are attached to a material point and are
determined by a local cartesian coordinates system (`mcsys`).

The integral is with respect to the volume of the domain ``V`` (i.e. a three
dimensional integral).

# Arguments
- `self` = finite element machine;
- `assembler` = assembler of the global matrix;
- `geom` = geometry field;
- `u` = nodal field to define the degree of freedom numbers;
- `cf`= data cache, which is called to evaluate the coefficient ``c``, given the
  location of the integration point, the Jacobian matrix, and the finite
  element label.
"""
function bilform_diffusion(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    u::NodalField{T},
    cf::DC,
) where {FEMM<:AbstractFEMM,A<:AbstractSysmatAssembler,FT,T,DC<:DataCache}
    if isempty(size(cf)) # scalar
        return _bilform_diffusion_iso(self, assembler, geom, u, cf)
    else
        return _bilform_diffusion_general(self, assembler, geom, u, cf)
    end
end

function _bilform_diffusion_general(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    u::NodalField{T},
    cf::DC,
) where {FEMM<:AbstractFEMM,A<:AbstractSysmatAssembler,FT,T,DC<:DataCache}
    fes = finite_elements(self)
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buff_b(self, geom, u)
    elmdim, elmat, elvec = _buff_e(self, geom, u, assembler)
    RmTJ, c_gradNT = _buff_d(self, geom, u)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    startassembly!(assembler, prod(size(elmat)) * count(fes), nalldofs(u), nalldofs(u))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elmat, zero(T)) # Initialize element matrix
        for j  in  1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            updatecsmat!(self.mcsys, loc, J, i, j)
            mulCAtB!(RmTJ, csmat(self.mcsys), J) # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], RmTJ)
            c = cf(loc, J, i, j)
            add_gkgt_ut_only!(elmat, gradN, (Jac * w[j]), c, c_gradNT)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(u, dofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function _bilform_diffusion_iso(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    u::NodalField{T},
    cf::DC,
) where {FEMM<:AbstractFEMM,A<:AbstractSysmatAssembler,FT,T,DC<:DataCache}
    fes = finite_elements(self)
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buff_b(self, geom, u)
    elmdim, elmat, elvec = _buff_e(self, geom, u, assembler)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    startassembly!(assembler, prod(size(elmat)) * count(fes), nalldofs(u), nalldofs(u))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elmat, zero(T)) # Initialize element matrix
        for j  in  1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            gradN!(fes, gradN, gradNparams[j], J)
            c = cf(loc, J, i, j)
            add_mggt_ut_only!(elmat, gradN, (c * Jac * w[j]))
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(u, dofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function bilform_diffusion(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    f::DC,
) where {FEMM<:AbstractFEMM,FT,T,DC<:DataCache}
    assembler = SysmatAssemblerSparseSymm()
    return bilform_diffusion(self, assembler, geom, u, f)
end

"""
    bilform_convection(
        self::FEMM,
        assembler::A,
        geom::NodalField{FT},
        u::NodalField{T},
        Q::NodalField{QT},
        rhof::DC
    ) where {FEMM<:AbstractFEMM, A<:AbstractSysmatAssembler, FT, T, QT, DC<:DataCache}

Compute the sparse matrix implied by the bilinear form of the "convection" type.

```math
\\int_{V}  {w} \\rho \\mathbf{u} \\cdot \\nabla q \\; \\mathrm{d} V
```

Here ``w`` is the scalar test function, ``\\mathbf{u}`` is the convective
velocity, ``q`` is the scalar trial function, ``\\rho`` is the mass density;
``\\rho`` is computed by `rhof`, which is a given function(data). Both test and
trial functions are assumed to be from the same approximation space. `rhof` is
represented with [`DataCache`](@ref), and needs to return a  scalar mass
density.

The integral is with respect to the volume of the domain ``V`` (i.e. a three
dimensional integral).

# Arguments
- `self` = finite element machine;
- `assembler` = assembler of the global matrix;
- `geom` = geometry field;
- `u` = convective velocity field;
- `Q` = nodal field to define the degree of freedom numbers;
- `rhof`= data cache, which is called to evaluate the coefficient ``\\rho``,
  given the location of the integration point, the Jacobian matrix, and the
  finite element label.
"""
function bilform_convection(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    u::NodalField{T},
    Q::NodalField{QT},
    rhof::DC,
) where {FEMM<:AbstractFEMM,A<:AbstractSysmatAssembler,FT,T,QT,DC<:DataCache}
    fes = finite_elements(self)
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buff_b(self, geom, Q)
    nsd = ndofs(u)
    eus = deepcopy(ecoords)
    elmdim, elmat, _ = _buff_e(self, geom, Q, assembler)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    startassembly!(assembler, prod(size(elmat)) * count(fes), nalldofs(Q), nalldofs(Q))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        gathervalues_asmat!(u, eus, fes.conn[i])
        fill!(elmat, zero(T)) # Initialize element matrix
        for j  in  1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            gradN!(fes, gradN, gradNparams[j], J)
            rho = rhof(loc, J, i, j)
            for p = 1:nne
                for r = 1:nne
                    accum = zero(eltype(elmat))
                    for s = 1:nsd
                        u_s = zero(T)
                        for q = 1:nne
                            u_s += Ns[j][q] * eus[q, s]
                        end
                        accum += u_s * gradN[r, s]
                    end
                    elmat[p, r] += Ns[j][p] * accum * (Jac * w[j])
                end
            end
        end # Loop over quadrature points
        gatherdofnums!(Q, dofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)# assemble matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function bilform_convection(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    Q::NodalField{QT},
    rhof::DC,
) where {FEMM<:AbstractFEMM,FT,T,QT,DC<:DataCache}
    assembler = SysmatAssemblerSparse() # must be able to assem. unsymmetric mtx
    return bilform_convection(self, assembler, geom, u, Q, rhof)
end

"""
    bilform_div_grad(
        self::FEMM,
        assembler::A,
        geom::NodalField{FT},
        u::NodalField{T},
        viscf::DC
    ) where {FEMM<:AbstractFEMM, A<:AbstractSysmatAssembler, FT, T, DC<:DataCache}

Compute the sparse matrix implied by the bilinear form of the "div grad" type.

```math
\\int_{V}  \\mu \\nabla \\mathbf{w}:  \\nabla\\mathbf{u}   \\; \\mathrm{d} V
```

Here `` \\mathbf{w}`` is the vector test function, ``\\mathbf{u}`` is the
velocity, ``\\mu`` is the dynamic viscosity (or kinematic viscosity, depending
on the formulation); ``\\mu`` is computed by `viscf`, which is a given function
(data). Both test and trial functions are assumed to be from the same
approximation space. `viscf` is represented with [`DataCache`](@ref), and needs
to return a  scalar viscosity.

The integral is with respect to the volume of the domain ``V`` (i.e. a three
dimensional integral).

# Arguments
- `self` = finite element machine;
- `assembler` = assembler of the global matrix;
- `geom` = geometry field;
- `u` = velocity field;
- `viscf`= data cache, which is called to evaluate the coefficient ``\\mu``,
  given the location of the integration point, the Jacobian matrix, and the
  finite element label.
"""
function bilform_div_grad(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    u::NodalField{T},
    viscf::DC,
) where {FEMM<:AbstractFEMM,A<:AbstractSysmatAssembler,FT,T,DC<:DataCache}
    fes = finite_elements(self)
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buff_b(self, geom, u)
    elmdim, elmat, _ = _buff_e(self, geom, u, assembler)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    startassembly!(assembler, prod(size(elmat)) * count(fes), nalldofs(u), nalldofs(u))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elmat, zero(T)) # Initialize element matrix
        for j  in  1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            gradN!(fes, gradN, gradNparams[j], J)
            mu = viscf(loc, J, i, j)
            factor = mu * (Jac * w[j])
            for a = 1:nne
                for b = 1:nne
                    for s = 1:ndn
                        p = (a - 1) * ndn + s
                        r = (b - 1) * ndn + s
                        for q = 1:ndn
                            elmat[p, r] += factor * gradN[a, q] * gradN[b, q]
                        end
                        for q = 1:ndn
                            r = (b - 1) * ndn + q
                            elmat[p, r] += factor * gradN[a, q] * gradN[b, s]
                        end
                    end
                end
            end
        end # Loop over quadrature points
        gatherdofnums!(u, dofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)# assemble matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function bilform_div_grad(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    viscf::DC,
) where {FEMM<:AbstractFEMM,FT,T,DC<:DataCache}
    assembler = SysmatAssemblerSparse()
    return bilform_div_grad(self, assembler, geom, u, viscf)
end

function _buff_cb(self, geom, u, mr, assembler)
    fes = finite_elements(self)
    ndn = ndofs(u) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes per element
    mdim = manifdim(fes) # manifold dimension of the element
    FT = eltype(geom.values)
    mdim = manifdim(fes) # manifold dimension of the element
    nstrs = nstressstrain(mr)  # number of stresses
    elmatdim = ndn * nne             # dimension of the element matrix
    B = fill(zero(FT), nstrs, elmatdim) # strain-displacement matrix -- buffer
    CB = fill(zero(FT), nstrs, elmatdim) # strain-displacement matrix -- buffer
    RmTJ = fill(zero(FT), mdim, mdim) # buffer
    return B, CB, RmTJ
end

"""
    bilform_lin_elastic(
        self::FEMM,
        assembler::A,
        geom::NodalField{FT},
        u::NodalField{T},
        cf::DC
    ) where {FEMM<:AbstractFEMM, A<:AbstractSysmatAssembler, FT, T, DC<:DataCache}

Compute the sparse matrix implied by the bilinear form of the "linearized elasticity" type.

```math
\\int_{V} (B \\mathbf{w})^T C  B \\mathbf{u}   \\; \\mathrm{d} V
```

Here `` \\mathbf{w}`` is the vector test function, ``\\mathbf{u}`` is the
displacement (velocity), ``C`` is the elasticity (viscosity) matrix; ``C`` is
computed by `cf`, which is a given function(data). Both test and trial
functions are assumed to be from the same approximation space. `cf` is
represented with [`DataCache`](@ref), and needs to return a matrix of the
appropriate size.

The integral is with respect to the volume of the domain ``V`` (i.e. a three
dimensional integral).

# Arguments
- `self` = finite element machine;
- `assembler` = assembler of the global matrix;
- `geom` = geometry field;
- `u` = velocity field;
- `viscf`= data cache, which is called to evaluate the coefficient ``\\mu``,
  given the location of the integration point, the Jacobian matrix, and the
  finite element label.
"""
function bilform_lin_elastic(
    self::FEMM,
    assembler::A,
    geom::NodalField{FT},
    u::NodalField{T},
    mr::Type{MR},
    cf::DC,
) where {
    FEMM<:AbstractFEMM,
    A<:AbstractSysmatAssembler,
    FT,
    T,
    MR<:AbstractDeforModelRed,
    DC<:DataCache,
}
    fes = finite_elements(self)
    nne, ndn, ecoords, dofnums, loc, J, gradN = _buff_b(self, geom, u)
    elmdim, elmat, _ = _buff_e(self, geom, u, assembler)
    B, CB, RmTJ = _buff_cb(self, geom, u, mr, assembler)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    startassembly!(assembler, prod(size(elmat)) * count(fes), nalldofs(u), nalldofs(u))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(elmat, zero(T)) # Initialize element matrix
        for j  in  1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            updatecsmat!(self.mcsys, loc, J, i, j)
            At_mul_B!(RmTJ, csmat(self.mcsys), J) # local Jacobian matrix
            gradN!(fes, gradN, gradNparams[j], RmTJ)
            blmat!(mr, B, Ns[j], gradN, loc, csmat(self.mcsys))
            C = cf(loc, J, i, j)
            add_btdb_ut_only!(elmat, B, Jac * w[j], C, CB)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(u, dofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums)# assemble matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function bilform_lin_elastic(
    self::FEMM,
    geom::NodalField{FT},
    u::NodalField{T},
    mr::Type{MR},
    cf::DC,
) where {FEMM<:AbstractFEMM,FT,T,MR<:AbstractDeforModelRed,DC<:DataCache}
    assembler = SysmatAssemblerSparseSymm()
    return bilform_lin_elastic(self, assembler, geom, u, mr, cf)
end

end # module



