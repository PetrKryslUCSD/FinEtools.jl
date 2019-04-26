"""
    FENodeSetModule

Module for the finite element node set.
"""
module FENodeSetModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import Base.count

"""
Finite element node set type.
"""
mutable struct FENodeSet
	"""
	Array of node locations. Indexed with the node number.
	The location of node j is given by `xyz[j,:]`.
	"""
	xyz::Array{FFlt, 2}

	function FENodeSet(xyz::FFltMat)
		self = new(deepcopy(xyz)); # Need to make a COPY of the input array!
		return self
	end
end


"""
    spacedim(self::FENodeSet)

Number of dimensions of the space in which the node lives, 1, 2, or 3.
"""
spacedim(self::FENodeSet) =size(self.xyz,2)

"""
    xyz3(self::FENodeSet)

Get the  3-D coordinate that define the location  of the node.
Even if the nodes  were specified in  lower dimension (1-D, 2-D)
this function returns  a 3-D coordinate  by padding with zeros.
"""
function xyz3(self::FENodeSet)
    if (size(self.xyz,2)==1)
        val = [self.xyz zeros(size(self.xyz,1),2)];
    elseif (size(self.xyz,2)==2)
        val = [self.xyz zeros(size(self.xyz,1),1)];
    else
        val = deepcopy(self.xyz);
    end
end

"""
    count(self::FENodeSet)

Get the number of finite element nodes in the node set.
"""
function count(self::FENodeSet)
    return size(self.xyz,1)::FInt
end


end
