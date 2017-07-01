module FENodeToFEMapModule


using FinEtools.FTypesModule

  # Map from finite element nodes to the finite elements.
    #
    # Map from finite element nodes to the finite elements connecting them.
        # For each  node referenced in the connectivity of
        # the finite element set on input, the numbers of the individual
        # finite elements that reference that node is stored in an array in
        # the so array map{}.
        #         Example: fes.conn= [7,6,5;
        #                             4,1,3;
        #                             3,7,5];
        #             The map reads
        #                     map{1} = [2];
        #                     map{2} = [];#  note that node number 2 is not referenced by the connectivity
        #                     map{3} = [2,3];
        #                     map{4} = [2];
        #                     map{5} = [1,3];
        #                     map{6} = [1];
        #                     map{7} = [1,3];
        # The individual elements from the connectivity that reference
        # node number 5 are 1 and 3, so that fes.conn(map{5},:) lists all the
        # nodes that are connected to node 5 (including node 5 itself).
struct FENodeToFEMap
    map::Array{FIntVec,1}
end
export FENodeToFEMap


function FENodeToFEMap(conns::FIntMat,nmax::FInt)
    #n=  maximum(conns[:]);#largest node number
    map = Array{FIntVec}(nmax)
    for i = 1:length(map)
        map[i] = []  # initially empty arrays
    end

    for i = 1:size(conns,2)
        for j = 1:size(conns,1)
            ni = conns[j,i];
            push!(map[ni],j)
        end
    end
    return FENodeToFEMap(map)
end

end
