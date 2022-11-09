"""
    MeshImportModule

Module for mesh import  functions.
"""
module MeshImportModule

__precompile__(true)

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import ..FENodeSetModule: FENodeSet
import ..FESetModule: setlabel!, nodesperelem
import ..FESetModule: AbstractFESet, FESetL2, FESetL3, FESetP1
import ..FESetModule: FESetQ4, FESetQ8, FESetQ9
import ..FESetModule: FESetT3, FESetT6, FESetT4, FESetT10
import ..FESetModule: FESetH8, FESetH20, FESetH27
import ..MeshModificationModule: renumberconn!
import ..MeshExportModule: H5MESH
import ..MeshExportModule: VTK
import Unicode: uppercase, isdigit
import LinearAlgebra: norm
using DelimitedFiles

"""
!!! note

    The arrays are reallocated as the data is read from files. The size of the
    arrays is increased by this much. If the number of entities to be read  is
    large, the chunk should be probably increased so that only a few reallocations
    are needed.
"""
const chunk = 1000


# Fix up an old style floating-point number without the exponent letter.  
function _fixupdecimal(s)
    os = ""
    for i = length(s):-1:1
        os = s[i] * os
        if (s[i] == '-') && (i > 1) && (uppercase(s[i-1]) != "E") && isdigit(s[i-1])
            os = "E" * os
        end
        if (s[i] == '+') && (i > 1) && (uppercase(s[i-1]) != "E") && isdigit(s[i-1])
            os = "E" * os
        end
    end
    return os
end 

"""
    import_NASTRAN(filename)

Import tetrahedral (4- and 10-node) NASTRAN mesh (.nas file).

Limitations:
1. only the GRID and CTETRA  sections are read.
2. Only 4-node and 10-node tetrahedra  are handled.
3. The file should be free-form (data separated by commas). 
Some fixed-format files can also be processed (large-field, but not small-field).

# Output
Data dictionary, with keys 
- "`fens`" (finite element nodes), 
- "`fesets`" (array of finite element sets).
"""
function import_NASTRAN(filename; allocationchunk=chunk, expectfixedformat = false)
    lines = readlines(filename)

    nnode = 0
    node = zeros(allocationchunk, 4)
    maxnodel = 10
    nelem = 0
    elem = zeros(FInt, allocationchunk, maxnodel + 3)

    current_line = 1
    while true
        if current_line >= length(lines)
            break
        end
        temp  = lines[current_line]
        current_line = current_line + 1
        if (length(temp) >= 4) && (uppercase(temp[1:4]) == "GRID")
            fixedformat = expectfixedformat
            largefield = false
            if (length(temp) >= 5) && (uppercase(temp[1:5]) == "GRID*")
                largefield = true; fixedformat = true
            end
            @assert (!fixedformat) || (fixedformat && largefield) "Can handle either free format or large-field fixed format"
            nnode = nnode + 1
            if size(node, 1) < nnode
                node = vcat(node, zeros(allocationchunk, 4))
            end
            if fixedformat
                # $------1-------2-------3-------4-------5-------6-------7-------8-------9-------0
                # GRID*                  5               0-1.66618812195+1-3.85740337853+0
                # *       1.546691269367+1               0
                node[nnode, 1] = parse(Float64, _fixupdecimal(temp[9:24]))
                node[nnode, 2] = parse(Float64, _fixupdecimal(temp[41:56]))
                node[nnode, 3] = parse(Float64, _fixupdecimal(temp[57:72]))
                temp  = lines[current_line]
                current_line = current_line + 1
                node[nnode, 4] = parse(Float64, _fixupdecimal(temp[9:24]))
            else  
                # Template:
                #   GRID,1,,-1.32846E-017,3.25378E-033,0.216954
                A = split(replace(temp, "," => " "))
                for  six = 1:4
                    node[nnode, six] = parse(Float64, A[six+1])
                end
            end
        end # GRID
        if (length(temp) >= 6) && (uppercase(temp[1:6]) == "CTETRA")
            # Template:
            #   CTETRA,1,3,15,14,12,16,8971,4853,8972,4850,8973,4848
            nelem = nelem + 1
            if size(elem, 1) < nelem
                elem = vcat(elem, zeros(FInt, allocationchunk, maxnodel + 3))
            end
            continuation = ""
            fixedformat = (length(temp) == 72) || expectfixedformat # Is this fixed-format record? This is a guess based on the length of the line.
            if fixedformat
                continuation = temp[min(length(temp), 72):length(temp)]
                temp = temp[1:min(length(temp), 72)]
            end
            A = split(replace(temp, "," => " "))
            elem[nelem, 1] = parse(FInt, A[2])
            elem[nelem, 2] = parse(FInt, A[3])
            if length(A) == 7  #  nodes per element  equals  4
                nperel = 4
            else
                nperel = 10
                if length(A) < 13 # the line is continued: read the next line
                    temp  = lines[current_line]
                    current_line = current_line + 1
                    if fixedformat
                        temp = temp[length(continuation)+1:min(length(temp), 72)]
                    end
                    temp = strip(temp)
                    Acont = split(replace(temp, "," => " "))
                    A = vcat(A, Acont)
                end
            end
            elem[nelem, 3] = nperel
            for  six = 1:nperel
                elem[nelem, six+3] = parse(FInt, A[six+3])
            end
        end # CTETRA

    end # while

    # Trim the arrays
    node = node[1:nnode,:]
    elem = elem[1:nelem,:]

    # The nodes need to be in serial order:  if they are not,  the element
    # connectivities  will not point at the right nodes
    @assert  norm(collect(1:nnode)-node[:,1]) == 0 "Nodes are not in serial order"

    # Process output arguments
    # Extract coordinates
    xyz=node[:,2:4]
    # Cleanup element connectivities
    ennod = unique(elem[:,3])
    @assert length(ennod) == 1 "Cannot handle mixture of element types"

    @assert ((ennod[1] == 4) || (ennod[1] == 10)) "Unknown element type"
    conn = elem[:,4:3+convert(FInt, ennod[1])]
    pids = elem[:,2] # property identifier: different identifiers, different materials
    upids = unique(pids)

    # Create output arguments. First the nodes
    fens = FENodeSet(xyz)

    # Now for each distinct property identifier, create one finite element set
    fesets = AbstractFESet[]
    for pid in upids
        ix = [i for i in 1:length(pids) if pids[i] == pid] 
        C = conn[ix, :]
        if ennod[1] == 4
            fes = FESetT4(C)
        else
            fes = FESetT10(C)
        end
        push!(fesets, fes)
    end
    
    output = FDataDict("fens"=>fens, "fesets"=>fesets, "property_ids"=>upids)
    return output
end

mutable struct _AbaqusElementSection
    ElementLine::AbstractString
    nelem::FInt
    elem::Array{FInt,2}
end

mutable struct _AbaqusNSetSection
    NSetLine::AbstractString
    nsetname::AbstractString
    generate::Bool
    nodes::Vector{FInt}
end

mutable struct _AbaqusElSetSection
    ESetLine::AbstractString
    esetname::AbstractString
    generate::Bool
    elements::Vector{FInt}
end

"""
    import_ABAQUS(filename)

Import Abaqus mesh (.inp file).

Limitations:
1. Only the `*NODE` and `*ELEMENT`  sections are read
2. Only 4-node and 10-node tetrahedra, 8-node or 20-node  hexahedra, 4-node
   quadrilaterals, 3-node triangles are handled.

# Output
Data dictionary, with keys 
- "`fens`" = finite element nodes.
- "`fesets`" = array of finite element sets.
- "`nsets`" = dictionary of "node sets", the keys are the names of the sets.
"""
function import_ABAQUS(filename; allocationchunk=chunk)
    lines = readlines(filename)

    maxelnodes = 20

    warnings = String[]
    nnode = 0
    node = zeros(allocationchunk, 4)
    Reading_nodes = false
    next_line = 1
    while true
        if next_line > length(lines)
            break
        end
        temp = uppercase(strip(lines[next_line]))
        next_line = next_line + 1
        if (length(temp) >= 5) && (temp[1:5] == "*NODE")
            Reading_nodes = true
            nnode = 0
            node = zeros(allocationchunk, 4)
            temp = uppercase(strip(lines[next_line]))
            next_line = next_line + 1
        end
        if Reading_nodes
            if temp[1:1] == "*" # another section started
                Reading_nodes = false
                break
            end
            nnode = nnode + 1
            if size(node, 1) < nnode # if needed, allocate more space
                node = vcat(node, zeros(allocationchunk, 4))
            end
            A = split(replace(temp, "," => " "))
            for  six = 1:length(A)
                node[nnode, six] = parse(Float64, A[six])
            end
        end
    end # while

    nelemset = 0
    elemset = _AbaqusElementSection[]
    Reading_elements = false
    next_line = 1
    while true
        if next_line > length(lines)
            break
        end
        temp = uppercase(strip(lines[next_line]))
        next_line = next_line + 1
        if (length(temp) >= 8) && (temp[1:8] == "*ELEMENT")
            Reading_elements = true
            nelemset = nelemset + 1
            nelem = 0
            a = _AbaqusElementSection(temp, nelem, zeros(FInt, allocationchunk, maxelnodes+1))
            push!(elemset, a)
            temp = uppercase(strip(lines[next_line]))
            next_line = next_line + 1
        end
        if Reading_elements
            if temp[1:1] == "*" # another section started
                Reading_elements = false
                break
            end
            elemset[nelemset].nelem = elemset[nelemset].nelem + 1
            if size(elemset[nelemset].elem, 1) < elemset[nelemset].nelem
                elemset[nelemset].elem = vcat(elemset[nelemset].elem,
                zeros(FInt, allocationchunk, maxelnodes+1))
            end
            A = split(temp, ",")
            if (A[end] == "") # the present line is continued on the next one
                temp = uppercase(strip(lines[next_line]))
                next_line = next_line + 1
                Acont = split(temp, ",")
                A = vcat(A[1:end-1], Acont)
            end
            for ixxxx = 1:length(A)
                elemset[nelemset].elem[elemset[nelemset].nelem, ixxxx] = parse(FInt, A[ixxxx])
            end
        end
    end # while

    # The nodes may not be in serial order. We need to renumber the nodes and adjust the connectivities.
    nodenumbers = FInt.(node[1:nnode, 1])
    newnumbering = fill(0, maximum(nodenumbers))
    newnumbering[nodenumbers] = 1:length(nodenumbers) 

    # truncate the array to just the lines read
    node = node[1:nnode, :] 
    

    nsetsarr = _AbaqusNSetSection[]
    Reading_nsetsarr = false
    next_line = 1
    while true
        if next_line > length(lines)
            break
        end
        temp = uppercase(strip(lines[next_line]))
        next_line = next_line + 1
        if (length(temp) >= 5) && (temp[1:5] == "*NSET")
            Reading_nsetsarr = true
            # *NSET, NSET=THE-REGION, GENERATE
            l = split(temp, ",")
            sn = ""
            generate = false
            for t in l
                t = rstrip(lstrip(t))
                ts = split(t, "=")
                if ts[1] == "NSET"
                    sn = ts[2]
                end
                if ts[1] == "GENERATE"
                    generate = true
                end
            end
            a = _AbaqusNSetSection(temp, sn, generate, FInt[])
            push!(nsetsarr, a)
            temp = uppercase(strip(lines[next_line]))
            next_line = next_line + 1
        end
        if Reading_nsetsarr
            if temp[1:1] == "*" # another section started
                Reading_nsetsarr = false
            else
                A = split(temp, ",")
                if nsetsarr[end].generate
                    fn = parse(FInt, A[1])
                    ln = parse(FInt, A[2])
                    st = parse(FInt, A[3])
                    nsetsarr[end].nodes = collect(fn:st:ln)
                else
                    for ixxxx = 1:length(A)
                        if !isempty(A[ixxxx])
                            nn = parse(FInt, A[ixxxx])
                            push!(nsetsarr[end].nodes,  nn)
                        end
                    end
                end
            end
        end
    end # while

    elsetsarr = _AbaqusElSetSection[]
    Reading_elsetsarr = false
    next_line = 1
    while true
        if next_line > length(lines)
            break
        end
        temp = uppercase(strip(lines[next_line]))
        next_line = next_line + 1
        if (length(temp) >= 6) && (temp[1:6] == "*ELSET")
            Reading_elsetsarr = true
            # *ELSET, ELSET=THE-REGION, GENERATE
            l = split(temp, ",")
            sn = ""
            generate = false
            for t in l
                t = rstrip(lstrip(t))
                ts = split(t, "=")
                if ts[1] == "ELSET"
                    sn = ts[2]
                end
                if ts[1] == "GENERATE"
                    generate = true
                end
            end
            a = _AbaqusElSetSection(temp, sn, generate, FInt[])
            push!(elsetsarr, a)
            temp = uppercase(strip(lines[next_line]))
            next_line = next_line + 1
        end
        if Reading_elsetsarr
            if temp[1:1] == "*" # another section started
                Reading_elsetsarr = false
            else
                A = split(temp, ",")
                if elsetsarr[end].generate
                    fn = parse(FInt, A[1])
                    ln = parse(FInt, A[2])
                    st = parse(FInt, A[3])
                    elsetsarr[end].elements = collect(fn:st:ln)
                else
                    for ixxxx = 1:length(A)
                        if !isempty(A[ixxxx])
                            en = parse(FInt, A[ixxxx])
                            push!(elsetsarr[end].elements,  en)
                        end
                    end
                end
            end
        end
    end # while

    # Process output arguments
    # Nodes
    xyz = node[:,2:4]
    fens = FENodeSet(xyz)

    # Element sets
    fesets = AbstractFESet[]

    function feset_construct(elemset1)
        temp = uppercase(strip(elemset1.ElementLine))
        b  = split(temp, ",")
        for ixxx = 1:length(b)
            c = split(b[ixxx], "=")
            if (uppercase(strip(c[1])) == "TYPE") && (length(c) > 1)
                TYPE = uppercase(strip(c[2]))
                if (length(TYPE) >= 4) && (TYPE[1:4] == "C3D8")
                    return FESetH8(elemset1.elem[:, 2:9])
                elseif (length(TYPE) >= 5) && (TYPE[1:5] == "C3D20")
                    return FESetH20(elemset1.elem[:, 2:21])
                elseif ((length(TYPE) >= 4) && (TYPE[1:4] == "C3D4")) || 
                    ((length(TYPE) >= 5) && (TYPE[1:5] == "DC3D4"))
                    return FESetT4(elemset1.elem[:, 2:5])
                elseif (length(TYPE) >= 5) && (TYPE[1:5] == "C3D10")
                    return FESetT10(elemset1.elem[:, 2:11])
                elseif ((length(TYPE) >= 2) && (TYPE[1:2] == "S3")) ||
                    ((length(TYPE) >= 3) && (TYPE[1:3] == "S3R"))  || 
                    ((length(TYPE) >= 5) && (TYPE[1:5] == "STRI3"))  || 
                    ((length(TYPE) >= 4) && (TYPE[1:4] == "CPS3")) || 
                    ((length(TYPE) >= 5) && (TYPE[1:5] == "DC2D3"))
                    return FESetT3(elemset1.elem[:, 2:4])
                elseif ((length(TYPE) >= 5) && (TYPE[1:5] == "DC2D4")) ||
                    ((length(TYPE) >= 4) && (TYPE[1:4] == "CPS4")) ||
                    ((length(TYPE) >= 2) && (TYPE[1:2] == "S4")) ||
                    ((length(TYPE) >= 3) && (TYPE[1:3] == "S4R"))
                    return FESetQ4(elemset1.elem[:, 2:5])
                elseif (length(TYPE) >= 4) && (TYPE[1:4] == "CPS6") 
                    return FESetT6(elemset1.elem[:, 2:7])
                else
                    return nothing
                end
            end
        end
    end

    for ixxxx = 1:length(elemset)
        elemset[ixxxx].elem = elemset[ixxxx].elem[1:elemset[ixxxx].nelem, :]
        fes = feset_construct(elemset[ixxxx])
        if (fes == nothing)
            push!(warnings, "Don't know how to handle " * elemset[ixxxx].ElementLine)
        else
            fes = renumberconn!(fes, newnumbering)
            push!(fesets, fes)
        end
    end

    nsets = Dict()
    for ixxxx = 1:length(nsetsarr)
        newnodes = newnumbering[nsetsarr[ixxxx].nodes]
        nsets[nsetsarr[ixxxx].nsetname] = newnodes
    end

    elsets = Dict()
    for ixxxx = 1:length(elsetsarr)
        elsets[elsetsarr[ixxxx].esetname] = elsetsarr[ixxxx].elements
    end

    output = FDataDict("fens"=>fens, "fesets"=>fesets, "nsets"=>nsets, "elsets"=>elsets, "warnings"=>warnings)
    return output
end

"""
    import_MESH(filename)

Import mesh in the MESH format (.mesh, .xyz, .conn triplet of files).

# Output
Data dictionary, with keys 
- "`fens`" = finite element nodes.
- "`fesets`" = array of finite element sets.
"""
function import_MESH(filename)
    warnings = String[]
    meshfilebase, ext = splitext(filename)
    if ext == ""
        ext = ".mesh"
    end
    meshfile = meshfilebase * ext
    meshfiledir = dirname(meshfile)
    # Mesh file in the format of the FinEtools .mesh file
    datinfo = open(meshfile, "r") do file
        readdlm(file)
    end
    # The first line is the name of the file with the coordinates
    Xfile = isfile(datinfo[1]) ? datinfo[1] : joinpath(meshfiledir, datinfo[1])
    X = open(Xfile, "r") do file
        readdlm(file, ' ', Float64)
    end
    # The second line is the name of the element name tag.
    etype = datinfo[2]
    Cfile = isfile(datinfo[3]) ? datinfo[3] : joinpath(meshfiledir, datinfo[3])
    C = open(Cfile, "r") do file
        readdlm(file, ' ', Int64)
    end

    fens = FENodeSet(X)

    if     etype == "H8"
        fes = FESetH8(C)
    elseif etype == "H20"
        fes = FESetH20(C)
    elseif etype == "T4"
        fes = FESetT4(C)
    elseif etype == "T10"
        fes = FESetT10(C)
    elseif etype == "Q4"
        fes = FESetQ4(C)
    elseif etype == "Q8"
        fes = FESetQ8(C)
    elseif etype == "T3"
        fes = FESetT3(C)
    elseif etype == "T6"
        fes = FESetT6(C)
    else
        push!(warnings, "Don't know how to handle " * etype)
        fes = nothing
    end

    # Is an optional label file present?
    if length(datinfo) >= 4
        labelfile = isfile(datinfo[4]) ? datinfo[4] : joinpath(meshfiledir, datinfo[4])
        lab = open(labelfile, "r") do file
            readdlm(file, ' ', Int64)
        end
        setlabel!(fes, vec(lab))
    end
  
    output = FDataDict("fens"=>fens, "fesets"=>[fes], "warnings"=>warnings)
    return output
end

using DataDrop

"""
    import_H5MESH(meshfile)

Import mesh in the H5MESH format (.h5mesh file).

# Output
Data dictionary, with keys 
- "`fens`" = finite element nodes.
- "`fesets`" = array of finite element sets.
"""
function import_H5MESH(meshfile)
    warnings = String[]
    meshfilebase, ext = splitext(meshfile)
    if ext == ""
        ext = ".h5mesh"
    end
    fname = DataDrop.with_extension(meshfile, ext) 
    
    X = DataDrop.retrieve_matrix(fname, "xyz")
    fens = FENodeSet(X)

    etype = DataDrop.retrieve_value(fname, "etype")

    C = DataDrop.retrieve_matrix(fname, "conn")
    if eltype(C) != FInt
        C = FInt.(C)
    end
    if     etype == VTK.H8
        fes = FESetH8(C)
    elseif etype == VTK.H20
        fes = FESetH20(C)
    elseif etype == VTK.T4
        fes = FESetT4(C)
    elseif etype == VTK.T10
        fes = FESetT10(C)
    elseif etype == VTK.Q4
        fes = FESetQ4(C)
    elseif etype == VTK.Q8
        fes = FESetQ8(C)
    elseif etype == VTK.T3
        fes = FESetT3(C)
    elseif etype == VTK.T6
        fes = FESetT6(C)
    else
        push!(warnings, "Don't know how to handle " * etype)
        fes = nothing
    end

    lab = DataDrop.retrieve_matrix(fname, "label")
    setlabel!(fes, vec(lab))
    
  
    output = FDataDict("fens"=>fens, "fesets"=>[fes], "warnings"=>warnings)
    return output
end

end #Module
