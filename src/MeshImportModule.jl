"""
    MeshImportModule

Module for mesh import  functions.  
"""
module MeshImportModule

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.FENodeSetModule

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

  chunk = 5000
  nnode = 0
  node = zeros(chunk, 4)
  nelem = 0
  elem = zeros(Int, chunk, 13)

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
        elem = vcat(elem, zeros(chunk, 4))
      end
      A = split(replace(temp, ",", " "))
      elem[nelem, 1] = parse(Int64, A[2])
      elem[nelem, 2] = parse(Int64, A[3])
      if length(A) == 7  #  nodes per element  equals  4
        nperel = 4
      else
        nperel = 10
        if length(A) < 13 # the line is continued: read the next line
          temp  = lines[current_line]
          current_line = current_line + 1
          temp = strip(temp)
          Acont = split(replace(temp, ",", " "))
          A = hcat(A, Acont)
        end
      end
      elem[nelem, 3] = nperel
      for  six = 1:nperel
        elem[nelem, six+3] = parse(Int64, A[six+3])
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
  conn = elem[:,4:3+ennod[1]]

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

end #Module
