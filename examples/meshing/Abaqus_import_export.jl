module junk

using FinEtools
using FinEtools.MeshExportModule

chunk = 50

mutable struct AbaqusElementSection
  ElementLine::AbstractString
  nelem::Int
  elem::Array{Int64}
end

"""
    import_ABAQUS(filename)

Import tetrahedral (4- and 10-node) Abaqus mesh (.inp file).

    Limitations:
    1. Only the `*NODE` and `*ELEMENT`  sections are read
    2. Only 4-node and 10-node tetrahedra  are handled.
    3. The file needs to be free-form (data separated by commas).

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
      for  six = 1:4
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
      a = AbaqusElementSection(temp, nelem, zeros(Int64, chunk, maxelnodes+1))
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
                                      zeros(Int64, chunk, maxelnodes+1))
      end
      A = split(temp, ",")
      if (A[end] == "") # the present line is continued on the next one
        temp = uppercase(strip(lines[next_line]))
        next_line = next_line + 1
        Acont = split(temp, ",")
        A = hcat(A[1:end-1], Acont)
      end
      for ixxxx = 1:length(A)
        elemset[nelemset].elem[elemset[nelemset].nelem, ixxxx] = parse(Int64, A[ixxxx])
      end
    end
  end # while

  node = node[1:nnode, :] # truncate the array to just the lines read
  # The nodes need to be in serial order:  if they are not,  the element
  # connectivities  will not point at the right nodes
  if norm(collect(1:nnode)-node[:,1]) != 0
    error("Nodes are not in serial order")
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
        else
          return nothing
        end
      end
    end
  end

  for ixxxx = 1:length(elemset)
    elemset[ixxxx].elem = elemset[ixxxx].elem[1:elemset[ixxxx].nelem, :]
    fes = feset_construct(elemset[ixxxx])
    !(fes === nothing) &&  push!(fesarray, fes)
  end

  return fens, fesarray
end


## Solid cylinder/taper/sphere—-temperature loading; quadratic brick mesh

# The mesh  will be created in a very coarse representation from the
# key points in the drawing. The first coordinate is radial, the second coordinate is axial.
rz=[1.     0.;#A
1.4    0.;#B
0.995184726672197   0.098017140329561;
1.393258617341076 0.137223996461385;
0.980785  0.195090;#
1.37309939 0.27312645;
0.956940335732209   0.290284677254462
1.339716470025092 0.406398548156247
0.9238795  0.38268;#C
1.2124  0.7;#D
0.7071  0.7071;#E
1.1062  1.045;#F
0.7071  (0.7071+1.79)/2;#(E+H)/2
1.      1.39;#G
0.7071  1.79;#H
1.      1.79;#I
]*phun("M")
tolerance =1.e-6*phun("M")

# This is the quadrilateral mesh of the cross-section.   It will be modified and
# refined as  we go.
fens = FENodeSet(rz);
fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);

nref = 1;
for ref = 1:nref
  fens,fes = Q4refine(fens,fes);
  list = selectnode(fens, distance=1.0+0.1/2^nref, from=[0. 0.], inflate=tolerance);
  fens.xyz[list,:] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
end

##
# The mesh is extruded by sweeping around the axis of symmetry.
# Only a single layer of elements is generated of internal angle
# |angslice|.
nLayers = 7;
angslice  = 5*pi/16;

##
# First the mesh is extruded to a block whose third dimension
# represents the angular coordinate.
fens,fes = H8extrudeQ4(fens, fes, nLayers,
(rz,k)->[rz[1],rz[2],0.0]-(k)/nLayers*[0.,0.,angslice]);
# The block is now converted  to the axially symmetric geometry by using the
# third (angular) coordinate  to sweep out  an axially symmetric domain. The
# ccoordinates of the nodes at this point are |rza|,  radial distance,
# Z-coordinate, angle.
sweep(rza) = [-rza[1]*sin(rza[3]+angslice/2.0), rza[1]*cos(rza[3]+angslice/2.0), rza[2]]
for j=1:size(fens.xyz,1)
  fens.xyz[j,:] = sweep(fens.xyz[j,:])
end

AE = AbaqusExporter("LE11NAFEMS_H8");
HEADING(AE, "LE11NAFEMS: Linear bricks.");
PART(AE, "part1");
END_PART(AE);
ASSEMBLY(AE, "ASSEM1");
INSTANCE(AE, "INSTNC1", "PART1");
NODE(AE, fens.xyz);
ELEMENT(AE, "c3d8rh", "AllElements", 1, fes.conn)
END_INSTANCE(AE);
END_ASSEMBLY(AE);
close(AE)


fens, fesarray = import_ABAQUS("./LE11NAFEMS_H8.inp")

File = "LE11NAFEMS_H8.vtk"
MeshExportModule.vtkexportmesh(File, fens, fesarray[1])
@async run(`"paraview.exe" $File`)



end
