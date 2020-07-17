 const GMSH_ELEMENTS = Dict(
 1=>FESetL2,
 2=>FESetT3, 
 3=>FESetQ4, 
 4=>FESetT4, 
 5=>FESetH8, 
 8=>FESetL3,
 9=>FESetT6,
 10=>FESetQ9,
 11=>FESetT10,
 12=>FESetH27,
 15=>FESetP1,
 16=>FESetQ8,
 17=>FESetH20,
 )

mutable struct GmshElementEntityBlock
    entityDim::Int64
    entityTag::Int64
    elementType::Int64
    numElementsInBlock::Int64
    elem::Array{FInt,2}
    nelem::Int64
    nnodes::Int64
end

"""
    import_ABAQUS(filename)

Import tetrahedral (4- and 10-node) or hexahedral (8- and 20-node) Gmsh mesh
(.inp file).

Limitations:
1. Only the `*NODE` and `*ELEMENT`  sections are read
2. Only 4-node and 10-node tetrahedra, 8-node or 20-node  hexahedra, 3-node triangles
    are handled.

# Output
Data dictionary, with keys 
- "`fens`" (finite element nodes), 
- "`fesets`" (array of finite element sets).
"""
function import_GMSH(filename; allocationchunk=chunk)
    lines = readlines(filename)

    warnings = String[]

    # $MeshFormat // same as MSH version 2
    #   version(ASCII double; currently 4.1)
    #     file-type(ASCII int; 0 for ASCII mode, 1 for binary mode)
    #     data-size(ASCII int; sizeof(size_t))
    #   < int with value one; only in binary mode, to detect endianness >
    # $EndMeshFormat
    version = 0.0
    file_type = 0
    data_size = 0

    next_line = 1
    while true
        if next_line > length(lines)
            break
        end
        temp = uppercase(strip(lines[next_line]))
        next_line = next_line + 1

        reading_meshformat = false
        if (length(temp) >= 11) && (temp[1:11] == "\$MESHFORMAT")
            reading_meshformat = true
            temp = uppercase(strip(lines[next_line]))
            next_line = next_line + 1
        end
        if reading_meshformat
            if temp[1:1] == "\$" # another section started
                reading_meshformat = false
                break
            end
            # temp = " 2.2  0  8"
            A = split(replace(temp, "," => " "))
            version = parse(Float64, A[1])
            if version < 3.0
                error("Cannot handle GMSH files of version < 3.0")
            end
            file_type = parse(Int64, A[2])
            if file_type != 0
                error("Cannot handle binary GMSH files")
            end
            data_size = parse(Int64, A[3])
        end
    end # while

    # $Nodes
    #   numEntityBlocks(size_t) numNodes(size_t)
    #     minNodeTag(size_t) maxNodeTag(size_t)
    #   entityDim(int) entityTag(int) parametric(int; 0 or 1)
    #     numNodesInBlock(size_t)
    #     nodeTag(size_t)
    #     ...
    #     x(double) y(double) z(double)
    #        < u(double; if parametric and entityDim >= 1) >
    #        < v(double; if parametric and entityDim >= 2) >
    #        < w(double; if parametric and entityDim == 3) >
    #     ...
    #   ...
    # $EndNodes
    # $Nodes
    # 1 6 1 6          1 entity bloc, 6 nodes total, min/max node tags: 1 and 6
    # 2 1 0 6          2D entity (surface) 1, no parametric coordinates, 6 nodes
    # 1                  node tag #1
    # 2                  node tag #2
    # 3                  etc.
    # 4
    # 5
    # 6
    # 0. 0. 0.           node #1 coordinates (0., 0., 0.)
    # 1. 0. 0.           node #2 coordinates (1., 0., 0.)
    # 1. 1. 0.           etc.
    # 0. 1. 0.
    # 2. 0. 0.
    # 2. 1. 0.
    # $EndNodes
    numEntityBlocks = 0
    numNodes = 0
    minNodeTag = 0
    maxNodeTag = 0
    nnode = 0
    node = zeros(allocationchunk, 4)
    reading_nodes = false
    reading_node_info = false
    reading_node_numbers = false
    reading_node_coordinates = false
    next_line = 1
    while true
        if next_line > length(lines)
            break
        end
        temp = uppercase(strip(lines[next_line]))
        next_line = next_line + 1

        if (length(temp) >= 6) && (temp[1:6] == "\$NODES")
            reading_node_info = true
            reading_nodes = true
            reading_node_numbers = false
            reading_node_coordinates = false
            nnode = 0
            node = zeros(allocationchunk, 4)
            temp = uppercase(strip(lines[next_line]))
            next_line = next_line + 1
        end
        if reading_nodes
            if reading_node_info 
                A = split(replace(temp, "," => " "))
                numEntityBlocks = parse(Int64, A[1])
                numNodes = parse(Int64, A[2])
                minNodeTag = parse(Int64, A[3])
                maxNodeTag = parse(Int64, A[4])
                if numEntityBlocks != 0
                    push!(warnings, "Don't know how to handle entities")
                end
                # Skip the entity info
                for s in 1:numEntityBlocks
                    next_line = next_line + 1
                end
                temp = uppercase(strip(lines[next_line]))
                next_line = next_line + 1
                reading_node_info = false
                reading_node_numbers = true
                reading_node_coordinates = false
                nnode = 0
            end
            if temp[1:1] == "\$" # another section started
                reading_nodes = false
                break
            end
            if reading_node_numbers
                nnode = nnode + 1
                if nnode > numNodes
                    reading_node_numbers = false
                    reading_node_coordinates = true
                    nnode = 0
                else
                    if size(node, 1) < nnode # if needed, allocate more space
                        node = vcat(node, zeros(allocationchunk, 4))
                    end
                    A = split(replace(temp, "," => " "))
                    six = 1
                    node[nnode, six] = parse(Float64, A[six])
                end
            end
            if reading_node_coordinates
                nnode = nnode + 1
                if nnode > numNodes
                    reading_node_numbers = false
                    reading_node_coordinates = false
                    nnode = 0
                else
                    A = split(replace(temp, "," => " "))
                    for  six in 1:min(length(A), 3)
                        # The first column are the node numbers
                        node[nnode, six+1] = parse(Float64, A[six])
                    end
                end
            end
            
        end
    end # while
    
    # $Elements
    #   numEntityBlocks(size_t) numElements(size_t)
    #     minElementTag(size_t) maxElementTag(size_t)
    #   entityDim(int) entityTag(int) elementType(int; see below)
    #     numElementsInBlock(size_t)
    #     elementTag(size_t) nodeTag(size_t) ...
    #     ...
    #   ...
    # $EndElements
    numEntityBlocks = 0
    numElements = 0
    minElementTag = 0
    maxElementTag = 0
    entityTag = 0
    elembl = GmshElementEntityBlock[]
    reading_elements = false
    reading_elements_ent_block = false
    reading_elements_data = false
    next_line = 1
    while true
        if next_line > length(lines)
            break
        end
        temp = uppercase(strip(lines[next_line]))
        next_line = next_line + 1
        if (length(temp) >= 9) && (temp[1:9] == "\$ELEMENTS")
            reading_elements = true
            temp = uppercase(strip(lines[next_line]))
            next_line = next_line + 1
            A = split(replace(temp, "," => " "))
            numEntityBlocks = parse(Int64, A[1])
            numElements = parse(Int64, A[2])
            minElementTag = parse(Int64, A[3])
            maxElementTag = parse(Int64, A[4])
            # @show numEntityBlocks, numElements, minElementTag, maxElementTag
            for i in 1:numEntityBlocks
                a = GmshElementEntityBlock(0, 0, 0, 0, zeros(FInt, 0, 0), 0, 0)
                push!(elembl, a)
            end
            # @show elembl
            temp = uppercase(strip(lines[next_line]))
            next_line = next_line + 1
            reading_elements_ent_block = true
        end
        if reading_elements
            @show next_line, temp  
            if temp[1:1] == "\$" # another section started
                reading_elements = false
                break
            end
            if reading_elements_ent_block
                A = split(temp, " ")
                entityDim = parse(Int64, A[1])
                entityTag = parse(Int64, A[2])
                elementType = parse(Int64, A[3])
                numElementsInBlock = parse(Int64, A[4])
                @show entityDim, entityTag, elementType, numElementsInBlock
                elembl[entityTag].entityTag = entityTag
                elembl[entityTag].entityDim = entityDim
                elembl[entityTag].elementType = elementType
                elembl[entityTag].numElementsInBlock = numElementsInBlock
                elembl[entityTag].nelem = 0
                elembl[entityTag].nnodes = nodesperelem(GMSH_ELEMENTS[elembl[entityTag].elementType])
                elembl[entityTag].elem = zeros(Int64, allocationchunk, elembl[entityTag].nnodes+1)
                temp = uppercase(strip(lines[next_line]))
                next_line = next_line + 1
                reading_elements_ent_block = false
                reading_elements_data = true
            end
            if reading_elements_data
                elembl[entityTag].nelem = elembl[entityTag].nelem + 1
                if size(elembl[entityTag].elem, 1) < elembl[entityTag].nelem
                    elembl[entityTag].elem = vcat(elembl[entityTag].elem,
                    zeros(FInt, allocationchunk, elembl[entityTag].nnodes+1))
                end
                A = split(temp, " ")
                for ixxxx = 1:(elembl[entityTag].nnodes+1)
                    elembl[entityTag].elem[elembl[entityTag].nelem, ixxxx] = parse(FInt, A[ixxxx])
                end
                # @show elembl[entityTag].nelem
                # @show elembl[entityTag].elem[elembl[entityTag].nelem, :]
                if elembl[entityTag].nelem == elembl[entityTag].numElementsInBlock
                    reading_elements_data = false
                    reading_elements_ent_block = entityTag < numEntityBlocks
                end
            end
        end
    end # while
    # @show elembl
@show node
    node = node[1:nnode, :] # truncate the array to just the lines read
    # The nodes need to be in serial order:  if they are not,  the element
    # connectivities  will not point at the right nodes. So,  if that's the case we
    # will figure out the sequential numbering of the nodes  and then we will
    # renumber the connectivities of the elements.
    newnumbering = collect(1:nnode)
    if norm(collect(1:nnode)-node[:,1]) != 0
        newnumbering = zeros(FInt, convert(FInt, maximum(node[:,1])))
        jn = 1
        for ixxxx = 1:size(node, 1)
            if node[ixxxx,1] != 0
                on = convert(FInt, node[ixxxx,1])
                newnumbering[on] = jn
                jn = jn + 1
            end
        end
    end

    # Process output arguments
    # Nodes
    xyz = node[:,2:4]
    fens = FENodeSet(xyz)

    # Element sets
    fesets = AbstractFESet[]

    function feset_construct(elembl1)
        typ = GMSH_ELEMENTS[elembl1.elementType]
        return typ(elembl1.elem[:, 2:end])
    end

    for ixxxx = 1:length(elembl)
        elembl[ixxxx].elem = elembl[ixxxx].elem[1:elembl[ixxxx].nelem, :]
        fes = feset_construct(elembl[ixxxx])
        if (fes == nothing)
            push!(warnings, "Don't know how to handle element type $(elembl[ixxxx].elementType)")
        else
            fes = renumberconn!(fes, newnumbering)
            push!(fesets, fes)
        end
    end

    output = FDataDict("fens"=>fens, "fesets"=>fesets, "warnings"=>warnings)
    return output
end
