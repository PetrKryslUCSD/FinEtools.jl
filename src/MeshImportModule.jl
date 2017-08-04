"""
    MeshImportModule

Module for mesh import  functions.
"""
module MeshImportModule

using FinEtools
using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.FENodeSetModule


const chunk = 50

"""
    import_NASTRAN(filename)

Import tetrahedral (4- and 10-node) NASTRAN mesh (.nas file).

    Limitations:
    1. only the GRID and CTETRA  sections are read.
    2. Only 4-node and 10-node tetrahedra  are handled.
    3. The file needs to be free-form (data separated by commas).

"""
function import_NASTRAN(filename)
  lines = readlines(filename)

  nnode = 0
  node = zeros(chunk, 4)
  maxnodel = 10
  nelem = 0
  elem = zeros(FInt, chunk, maxnodel + 3)

  current_line = 1
  while true
    if current_line >= length(lines)
      break
    end
    temp  = lines[current_line]
    current_line = current_line + 1
    temp = strip(temp)
    temp = uppercase(temp)
    if temp[1:4] == "GRID"
      # Template:
      #   GRID,1,,-1.32846E-017,3.25378E-033,0.216954
      nnode = nnode + 1
      if size(node, 1) < nnode
        node = vcat(node, zeros(chunk, 4))
      end
      A = split(replace(temp, ",", " "))
      for  six = 1:4
        node[nnode, six] = parse(Float64, A[six+1])
      end
    end # GRID
    if temp[1:6] == "CTETRA"
      # Template:
      #   CTETRA,1,3,15,14,12,16,8971,4853,8972,4850,8973,4848
      nelem = nelem + 1
      if size(elem, 1) < nelem
        elem = vcat(elem, zeros(FInt, chunk, maxnodel + 3))
      end
      A = split(replace(temp, ",", " "))
      elem[nelem, 1] = parse(FInt, A[2])
      elem[nelem, 2] = parse(FInt, A[3])
      if length(A) == 7  #  nodes per element  equals  4
        nperel = 4
      else
        nperel = 10
        if length(A) < 13 # the line is continued: read the next line
          temp  = lines[current_line]
          current_line = current_line + 1
          temp = strip(temp)
          Acont = split(replace(temp, ",", " "))
          A = vcat(A, Acont)
        end
      end
      elem[nelem, 3] = nperel
      for  six = 1:nperel
        elem[nelem, six+3] = parse(FInt, A[six+3])
      end
    end # CTETRA

  end # while

  node=node[1:nnode,:]
  elem=elem[1:nelem,:]
  # The nodes need to be in serial order:  if they are not,  the element
  # connectivities  will not point at the right nodes
  if norm(collect(1:nnode)-node[:,1]) != 0
    error("Nodes are not in serial order")
  end

    # Process output arguments
  # Extract coordinates
  xyz=node[:,2:4]
  # Cleanup element connectivities
  ennod = unique(elem[:,3])
  if  length(ennod) != 1
    error("This function cannot treat a mixture of element types at this point")
  end

  if (ennod[1] != 4) && (ennod[1] != 10)
    error("Unknown element type")
  end
  conn = elem[:,4:3+convert(FInt, ennod[1])]

  # Create output arguments. First the nodes
  fens = FENodeSet(xyz)
  if ennod[1] == 4
    fes = FESetT4(conn)
  else
    fes = FESetT10(conn)
  end
  setlabel!(fes, elem[:,2])

  return fens, fes
end

mutable struct AbaqusElementSection
  ElementLine::AbstractString
  nelem::FInt
  elem::Array{FInt,2}
end

"""
    import_ABAQUS(filename)

Import tetrahedral (4- and 10-node) or hexahedral (8- and 20-node) Abaqus mesh
(.inp file).

    Limitations:
    1. Only the `*NODE` and `*ELEMENT`  sections are read
    2. Only 4-node and 10-node tetrahedra and 8-node or 20-node  hexahedra
      are handled.

"""
function import_ABAQUS(filename)
  lines = readlines(filename)

  maxelnodes = 20

  nnode = 0
  node = zeros(chunk, 4)
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
      node = zeros(chunk, 4)
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
        node = vcat(node, zeros(chunk, 4))
      end
      A = split(replace(temp, ",", " "))
      for  six = 1:length(A)
        node[nnode, six] = parse(Float64, A[six])
      end
    end
  end # while

  nelemset = 0
  elemset = AbaqusElementSection[]
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
      a = AbaqusElementSection(temp, nelem, zeros(FInt, chunk, maxelnodes+1))
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
                                      zeros(FInt, chunk, maxelnodes+1))
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

  node = node[1:nnode, :] # truncate the array to just the lines read
  # The nodes need to be in serial order:  if they are not,  the element
  # connectivities  will not point at the right nodes. So,  if that's the case we
  # will figure out the sequential numbering of the nodes  and then we will
  # renumbered the connectivvities of the elements.
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
  fesarray = FESet[]

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
        elseif (length(TYPE) >= 4) && (TYPE[1:4] == "C3D4")
          return FESetT4(elemset1.elem[:, 2:5])
        elseif (length(TYPE) >= 5) && (TYPE[1:5] == "C3D10")
          return FESetT10(elemset1.elem[:, 2:11])
        else
          return nothing
        end
      end
    end
  end

  for ixxxx = 1:length(elemset)
    elemset[ixxxx].elem = elemset[ixxxx].elem[1:elemset[ixxxx].nelem, :]
    fes = feset_construct(elemset[ixxxx])
      if (fes === nothing)
        warn("Don't know how to handle " * elemset[ixxxx].ElementLine)
      else
        fes = renumberconn!(fes, newnumbering)
        push!(fesarray, fes)
      end
    end

  return fens, fesarray
end

end #Module
