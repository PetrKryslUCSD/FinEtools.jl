"""
    VoxelTetMeshingModule

Module for meshing of voxel data sets with tetrahedra.  
"""
module VoxelTetMeshingModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.VoxelBoxModule: VoxelBoxVolume, voxeldims
import FinEtools.MeshTetrahedronModule: T4voximg, tetv
import FinEtools.FESetModule: connasarray
import FinEtools.MeshModificationModule: interior2boundary, vertexneighbors, smoothertaubin
import FinEtools.TetRemeshingModule: coarsen
import LinearAlgebra: norm
import Statistics: mean

mutable struct  ElementSizeWeightFunction
    influenceweight::FFlt
    center::Vector{FFlt}
    influenceradius::FFlt
end

function weight(self::ElementSizeWeightFunction, xyz::Vector{FFlt})
    r = norm(xyz - self.center);
    return 1.0+(self.influenceweight)*(1.0-1.0/(1+exp(-3*(4.0*r/self.influenceradius-3))));
end
  
"""
    ImageMesher{CoordT,DataT}

Tetrahedral image mesher.  
"""
mutable struct  ImageMesher{CoordT,DataT}
    box::VoxelBoxVolume{CoordT,DataT}
    currentelementsize::FFlt
    emptyvoxel::DataT
    notemptyvoxel::Vector{DataT}
    elementsizeweightfunctions::Vector{ElementSizeWeightFunction}
    havemesh::Bool
    t::Array{Int, 2}
    v::Array{FFlt, 2}
    tmid::Vector{Int}
end

function ImageMesher(box::VoxelBoxVolume{CoordT,DataT}, emptyvoxel::DataT, notemptyvoxel::Vector{DataT}) where {CoordT,DataT} 
    return ImageMesher(box, 0.0, emptyvoxel, notemptyvoxel, ElementSizeWeightFunction[], false, fill(0, 0, 0), fill(zero(FFlt), 0, 0), Int[])
end

function evaluateweights(self::ImageMesher)
    vertex_weight = ones(FFlt, size(self.v,1));
    if !isempty(self.elementsizeweightfunctions)
        for i = 1:size(self.v,1)
            vertex_weight[i] = 1.0
            for jewf=1:length(self.elementsizeweightfunctions)
                vertex_weight[i] = max(vertex_weight[i], weight(self.elementsizeweightfunctions[jewf], vec(self.v[i, :])))
            end
        end
    end 
    return vertex_weight;
end

function coarsenvolume!(self::ImageMesher, desired_element_size::FFlt)
    vertex_weight = evaluateweights(self)
    self.t, self.v, self.tmid = coarsen(self.t, self.v, self.tmid; surface_coarsening = false, desired_ts = desired_element_size, vertex_weight = vertex_weight);
    return self;
end
  
function coarsensurface!(self::ImageMesher, desired_element_size::FFlt)
    vertex_weight = evaluateweights(self)
    self.t, self.v, self.tmid = coarsen(self.t, self.v, self.tmid; surface_coarsening = true, desired_ts = desired_element_size, vertex_weight = vertex_weight);
    return self;
end
 
function volumes!(V::Array{FFlt, 1}, v::Array{FFlt, 2}, t::Array{Int, 2})
    for i = 1:size(t, 1)
        v11, v12, v13 = v[t[i, 1], :]
        v21, v22, v23 = v[t[i, 2], :]
        v31, v32, v33 = v[t[i, 3], :]
        v41, v42, v43 = v[t[i, 4], :]
        V[i] = tetv(v11, v12, v13, v21, v22, v23, v31, v32, v33, v41, v42, v43) 
    end
    return V
end

function allnonnegativevolumes(v, t)
    for i = 1:size(t, 1)
        v11, v12, v13 = v[t[i, 1], :]
        v21, v22, v23 = v[t[i, 2], :]
        v31, v32, v33 = v[t[i, 3], :]
        v41, v42, v43 = v[t[i, 4], :]
        if tetv(v11, v12, v13, v21, v22, v23, v31, v32, v33, v41, v42, v43) < 0.0
             return false
        end 
    end
    return true
end

function smooth!(self::ImageMesher, npass::Int = 5)
    if !allnonnegativevolumes(self.v, self.t) # This test is not strictly necessary,  just  while the remeshing procedure is in flux
        error("shouldn't be here")
    end 
    V = zeros(size(self.t, 1)) # tetrahedron volumes

    # Find boundary vertices
    bv = falses(size(self.v, 1));
    f = interior2boundary(self.t, [1 3 2; 1 2 4; 2 3 4; 1 4 3])
    
    # find neighbors for the SURFACE vertices
    fvn = vertexneighbors(f, size(self.v, 1));
    # Smoothing considering only surface connections
    fv = FFltMat[]
    trialfv = deepcopy(self.v)
    for pass = 1:npass
        fv =  smoothertaubin(trialfv, fvn, bv, 1, 0.5, -0.5);
        V = volumes!(V, fv, self.t)
        for i = 1:length(V)
            if V[i] < 0.0 # smoothing is only allowed if it results in non-negative volume
                c = self.t[i, :]
                fv[c, :] = trialfv[c, :]# undo the smoothing
            end 
        end
        copyto!(trialfv, fv)
    end
    
    # find neighbors for the VOLUME vertices
    vn =  vertexneighbors(self.t, size(self.v, 1));
    # Smoothing considering all connections through the volume
    bv[vec(f)] .= true;
    v = FFltMat[]
    trialv = deepcopy(fv)
    for pass = 1:npass
        v = smoothertaubin(trialv, vn, bv, 1, 0.5, -0.5); 
        anynegative = true; chk=1
        while anynegative
            # println("Checking volumes $chk")
            V = volumes!(V, v, self.t)
            anynegative = false
            for i = 1:length(V)
                if V[i] < 0.0 # smoothing is only allowed if it results in non-negative volume
                    c = self.t[i, :]
                    v[c, :] = self.v[c, :]# undo the smoothing
                    anynegative = true
                end 
            end
            chk = chk+1
        end
        copyto!(trialv, v)
    end

    # # Correction of the vertices of the surface
    # copyto!(trialv, v)
    # for i= 1:length(fvn)
    #     if (!isempty(fvn[i]))
    #         v[i,:] = fv[i,:];
    #     end
    # end
    # V = volumes!(V, v, self.t)
    # for i = 1:length(V)
    #     if V[i] < 0.0 # smoothing is only allowed if it results in non-negative volume
    #         c = self.t[i, :]
    #         v[c, :] = trialv[c, :]# undo the smoothing
    #     end 
    # end
    if !allnonnegativevolumes(v, self.t) # This test is not strictly necessary,  just  while the remeshing procedure is in flux
        error("shouldn't be here")
    end 

    copyto!(self.v, v) # save the final result
    return self
end

"""
    mesh!(self::ImageMesher, stretch::FFlt = 1.2)

Perform a meshing step.  

If  no mesh exists,  the initial mesh is created; otherwise a coarsening
sequence of coarsen surface -> smooth -> coarsen volume -> smooth is performed.

After meshing the vertices, tetrahedra, and material identifiers,  can be retrieved
as `self.v`, `self.t`, and `self.tmid`.
"""
function mesh!(self::ImageMesher, stretch::FFlt = 1.2)
    if !self.havemesh
        fens, fes = T4voximg(self.box.data, vec([voxeldims(self.box)...]), self.notemptyvoxel)
        self.v = deepcopy(fens.xyz)
        self.t = connasarray(fes)
        self.tmid = deepcopy(fes.label)
        self.currentelementsize = mean(vec([voxeldims(self.box)...]))
        smooth!(self);
        self.havemesh = true
    else
        coarsensurface!(self, sqrt(1.0)*self.currentelementsize)
        smooth!(self);
        coarsenvolume!(self, sqrt(2.0)*self.currentelementsize);
        smooth!(self);
        self.currentelementsize = stretch * self.currentelementsize
    end 
    return self
end

"""
    volumes(self::ImageMesher)

Compute tetrahedral volumes in the current mesh.  
"""
function volumes(self::ImageMesher)
    V = zeros(size(self.t, 1))
    return volumes!(V, self.v, self.t)
end

end