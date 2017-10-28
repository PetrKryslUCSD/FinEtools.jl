module VoxelTetMeshingModule

using FinEtools.FTypesModule
using FinEtools.TetRemeshingModule
using FinEtools.VoxelBoxModule
using FinEtools.MeshTetrahedronModule
using FinEtools.MeshModificationModule

export ElementSizeWeightFunction, ImageMesher, mesh!

mutable struct  ElementSizeWeightFunction
    influenceweight::FFlt
    center::Vector{FFlt}
    influenceradius::FFlt
end

function weight(self::ElementSizeWeightFunction, xyz::Vector{FFlt})
    r = norm(xyz - self.center);
    return 1.0+(self.influenceweight)*(1.0-1.0/(1+exp(-3*(4.0*r/self.influenceradius-3))));
end
  
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
    return ImageMesher(box, 0.0, emptyvoxel, notemptyvoxel, Vector{ElementSizeWeightFunction}(0), false, Array{Int, 2}(0, 0), Array{FFlt, 2}(0, 0), Vector{Int}(0))
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
    self.t, self.v, self.tmid = TetRemeshingModule.coarsen(self.t, self.v, self.tmid; surface_coarsening = false, desired_ts = desired_element_size, vertex_weight = vertex_weight);
    return self;
end
  
function coarsensurface!(self::ImageMesher, desired_element_size::FFlt)
    vertex_weight = evaluateweights(self)
    self.t, self.v, self.tmid = TetRemeshingModule.coarsen(self.t, self.v, self.tmid; surface_coarsening = true, desired_ts = desired_element_size, vertex_weight = vertex_weight);
    return self;
end
 
function smooth!(self::ImageMesher, npass::Int = 5)
    # fetmids = unique(self.tmid);
    # for j=1:length(fetmids)
    #     iix = find(m -> m==fetmids[j], self.tmid);
    #     self.v = MeshModificationModule.vsmoothing(self.v, self.t[iix,:]; npass = 5)
    # end 

    # Find boundary vertices
    bv = falses(size(self.v, 1));
    f = interior2boundary(self.t, [1 3 2; 1 2 4; 2 3 4; 1 4 3])
    
    # find neighbors for the SURFACE vertices
    fvn = MeshModificationModule.vertexneighbors(f, size(self.v, 1));
    # Smoothing considering only surface connections
    fv =  MeshModificationModule.smoothertaubin(self.v, fvn, bv, npass, 0.5, -0.5);
    
    # find neighbors for the VOLUME vertices
    vn =  MeshModificationModule.vertexneighbors(self.t, size(self.v, 1));
    # Smoothing considering all connections through the volume
    bv[vec(f)]=true;
    v = MeshModificationModule.smoothertaubin(self.v, vn, bv, npass, 0.5, -0.5); 

    # Correction of the vertices of the surface
    for i= 1:length(fvn)
        if (!isempty(fvn[i]))
            v[i,:] = fv[i,:];
        end
    end

    copy!(self.v, v)
    return self
end


function mesh!(self::ImageMesher, stretch::FFlt = 1.2)
    if !self.havemesh
        fens, fes = T4voximg(self.box.data, vec([voxeldims(self.box)...]), self.notemptyvoxel)
        self.v = deepcopy(fens.xyz)
        self.t = deepcopy(fes.conn)
        self.tmid = deepcopy(fes.label)
        self.currentelementsize = mean(vec([voxeldims(self.box)...]))
        smooth!(self);
        self.havemesh = true
    else
        coarsensurface!(self, sqrt(2.0)*self.currentelementsize)
        smooth!(self);
        coarsenvolume!(self, sqrt(3.0)*self.currentelementsize);
        smooth!(self);
        self.currentelementsize = stretch * self.currentelementsize
    end 
    return self
end

end