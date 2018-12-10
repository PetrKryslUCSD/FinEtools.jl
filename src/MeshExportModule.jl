"""
    MeshExportModule

Module for export of meshes and data defined on meshes.
"""
module MeshExportModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FESetModule: FESet, FESetP1, FESetL2, FESetT3, FESetQ4, FESetT4, FESetH8, FESetQ8, FESetL3, FESetT6, FESetT10, FESetH20, connasarray
import FinEtools.FENodeSetModule: FENodeSet
import Base.close

using Printf
import LinearAlgebra: norm, cross

################################################################################
# VTK export
################################################################################

const P1=1
const L2=3
const T3=5
const Q4=9
const T4=10
const H8=12
const L3=21
const T6=22
const Q8=23
const T10=24
const H20=25

VTKtypemap = Dict{DataType, Int}(FESetP1=>P1, FESetL2=>L2, FESetT3=>T3,
    FESetQ4=>Q4, FESetT4=>T4, FESetH8=>H8, FESetQ8=>Q8, FESetL3=>L3, FESetT6=>T6,
    FESetT10=>T10, FESetH20=>H20)

numnodesmap = Dict{Int, Int}(P1=>1, L2=>2, T3=>3,
    Q4=>4, T4=>4, H8=>8, Q8=>8, L3=>3, T6=>6,
    T10=>10, H20=>20)

"""
    vtkexportmesh(theFile::String, fens::FENodeSet, fes::T; opts...) where {T<:FESet}

Export mesh to a VTK 1.0 file as an unstructured grid.

`opts` = keyword argument list, where
`scalars` = array of tuples, (name, data)
`vectors` = array of tuples, (name, data)

For the `scalars`: If `data` is a vector, that data is exported as a single field.
On the other hand, if it is an 2d array, each column is exported  as a separate field.

"""
function vtkexportmesh(theFile::String, fens::FENodeSet, fes::T;
    opts...) where {T<:FESet}
    Cell_type = get(()->error("Unknown VTK type!"), VTKtypemap, typeof(fes));
    return vtkexportmesh(theFile, connasarray(fes), fens.xyz, Cell_type; opts...)
end

"""
    vtkexportmesh(theFile::String, Connectivity, Points, Cell_type;
        vectors=nothing, scalars=nothing)

Export mesh to a VTK 1.0 file as an unstructured grid.

`opts` = keyword argument list, where
`scalars` = array of tuples, (name, data)
`vectors` = array of tuples, (name, data)

For the `scalars`: If `data` is a vector, that data is exported as a single field.
On the other hand, if it is an 2d array, each column is exported  as a separate field.
"""
function vtkexportmesh(theFile::String, Connectivity, Points, Cell_type;
    vectors=nothing, scalars=nothing)
    X = Points
    if size(Points, 2) < 3
        X = zeros(size(Points,1),3)
        for j  = 1:size(Points,1)
            X[j,1:size(Points,2)] =  Points[j,:]
        end
    end
    numnodes = get(()->error("Wrong number of connected nodes!"), numnodesmap, Cell_type);
    if typeof(Connectivity[1]) <: Tuple # Vector of connectivity tuples: convert to an array
        c = fill(0, length(Connectivity), length(Connectivity[1]))
        for i = 1:length(Connectivity)
            c[i, :] = [Connectivity[i]...]
        end
        Connectivity = c
    end
    @assert numnodes == size(Connectivity,2)

    fid=open(theFile,"w");
    if (fid == -1)
        error(["Could not open " * theFile])
        return nothing
    end
    print(fid,"# vtk DataFile Version 1.0\n");
    print(fid,"FinEtools Export\n");
    print(fid,"ASCII\n");
    print(fid,"\n");
    print(fid,"DATASET UNSTRUCTURED_GRID\n");
    print(fid,"POINTS ", size(X,1), " double\n");
    for i= 1:size(X, 1)
        for j= 1:size(X,2)-1
            print(fid,X[i,j]," ");
        end
        print(fid,X[i,end],"\n");
    end
    print(fid,"\n");
    print(fid,"\n");

    print(fid,"CELLS ",size(Connectivity,1)," ",(size(Connectivity,1)*(size(Connectivity,2)+1)),"\n");
    for i= 1:size(Connectivity, 1)
        print(fid,size(Connectivity,2)," ");
        for j= 1:size(Connectivity,2)-1
            print(fid,Connectivity[i,j]-1," ");
        end
        print(fid,Connectivity[i,end]-1,"\n");
    end
    print(fid,"\n");
    print(fid,"\n");
    print(fid,"CELL_TYPES ",size(Connectivity,1),"\n");
    for i= 1:size(Connectivity,1)
        print(fid,Cell_type,"\n");
    end
    print(fid,"\n");
    print(fid,"\n");


    did_point_data = false

    # First try to write point data
    if (scalars != nothing)
        did_point_data = false
        for  sx = 1:length(scalars)
            name, data = scalars[sx]
            if (size(data, 1) == size(Points, 1)) # point data
                if (!did_point_data)
                    print(fid,"POINT_DATA ",size(data, 1),"\n");
                    did_point_data = true
                end
                if size(data, 2) > 1 # there are multiple scalar fields here
                    for i = 1:size(data, 2)
                        print(fid,"SCALARS ", name * "$(i)"," double\n");
                        print(fid,"LOOKUP_TABLE default\n");
                        for j= 1:size(data,  1)
                            print(fid,data[j, i],"\n");
                        end
                    end
                else  # there's just one scalar field here
                    print(fid,"SCALARS ", name," double\n");
                    print(fid,"LOOKUP_TABLE default\n");
                    for j= 1:size(data,  1)
                        print(fid,data[j],"\n");
                    end
                end
                print(fid,"\n");
            elseif (size(data, 1) == size(Connectivity, 1)) # cell data
                # Handled below, in the section for cell data
            end
        end
    end

    print(fid,"\n");

    if vectors != nothing
        for  vx = 1:length(vectors)
            name, data = vectors[vx]
            if (!did_point_data)
                print(fid,"POINT_DATA ",size(data, 1),"\n");
                did_point_data = true
            end
            print(fid,"VECTORS ", name," double\n");
            #print(fid,"LOOKUP_TABLE default\n");
            if size(data, 2) < 3
                X = zeros(size(data,1),3)
                for j  = 1:size(data,1)
                    X[j, 1:size(data,2)] =  data[j,:]
                end
            else
                X = data
            end
            for j= 1:size(X,1)
                k=1;
                print(fid,X[j,k]);
                for k = 2:size(X,2)
                    print(fid," ", X[j,k]);
                end
                print(fid,"\n");
            end
            print(fid,"\n");
        end
    end
    print(fid,"\n");

    # Are there any cell data?  If so, write out.
    if (scalars != nothing)
        did_cell_data = false
        for  sx = 1:length(scalars)
            name, data = scalars[sx]
            if (size(data, 1) == size(Points, 1)) # point data
                # Handled above for the case of point data
            elseif (size(data, 1) == size(Connectivity, 1)) # cell data
                if (!did_cell_data)
                    print(fid,"CELL_DATA ",size(data, 1),"\n");
                    did_cell_data = true
                end
                print(fid,"SCALARS ", name," double\n");
                print(fid,"LOOKUP_TABLE default\n");
                for j= 1:size(data,  1)
                    print(fid,data[j],"\n");
                end
                print(fid,"\n");
            end
        end
    end

    fid=close(fid);
    return true
end

"""
    vtkexportvectors(theFile::String, Points, vectors)

Export vector data to a VTK 1.0 file.

- `Points` = array of collection of coordinates (tuples or vectors), 
- `vectors` = array of tuples, `(name, data)`, where `name` is a string, and
  `data` is array of collection of coordinates (tuples or vectors).

# Example
```
Points = [(1.0, 3.0), (0.0, -1.0)]
vectors = [("v", [(-1.0, -2.0), (1.0, 1.0)])]
vtkexportvectors("theFile.VTK", Points, vectors)
```
will produce file with
```
# vtk DataFile Version 1.0
FinEtools Export
ASCII

DATASET UNSTRUCTURED_GRID
POINTS 2 double
1.0 3.0 0.0
0.0 -1.0 0.0


POINT_DATA 2
VECTORS v double
-1.0 -2.0 0.0
1.0 1.0 0.0
```
"""
function vtkexportvectors(theFile::String, Points, vectors)
    fid=open(theFile,"w");
    if (fid == -1)
        error(["Could not open " * theFile])
        return nothing
    end
    print(fid,"# vtk DataFile Version 1.0\n");
    print(fid,"FinEtools Export\n");
    print(fid,"ASCII\n");
    print(fid,"\n");
    print(fid,"DATASET UNSTRUCTURED_GRID\n");
    print(fid,"POINTS ", length(Points), " double\n");
    pt = fill(0.0, 3)
    for i = 1:length(Points)
        for k = 1:length(Points[i])
            pt[k] = Points[i][k]
        end
        for j= 1:length(pt)-1
            print(fid, pt[j]," ");
        end
        print(fid,pt[end],"\n");
    end
    print(fid,"\n");
    print(fid,"\n");

    did_point_data = false

    if vectors != nothing
        for  vx = 1:length(vectors)
            name, data = vectors[vx]
            if (!did_point_data)
                print(fid,"POINT_DATA ",length(data),"\n");
                did_point_data = true
            end
            print(fid,"VECTORS ", name," double\n");
            #print(fid,"LOOKUP_TABLE default\n");
            for i= 1:length(data)
                for k = 1:length(data[i])
                    pt[k] = data[i][k]
                end
                for j= 1:length(pt)-1
                   print(fid, pt[j]," ");
                end
                print(fid,pt[end],"\n");
            end
            print(fid,"\n");
        end
    end
    print(fid,"\n");

    fid=close(fid);
    return true
end

################################################################################
# Abaqus export
################################################################################

"""
    AbaqusExporter

Export mesh to Abaqus.
"""
mutable struct AbaqusExporter
  filename::AbstractString
  ios::IO
  element_range::Tuple{Int64, Int64}
  function AbaqusExporter(filename::AbstractString)
    if match(r".*\.inp$", filename) == nothing
      filename = filename * ".inp"
    end
    ios = open(filename,"w+")
    return new(deepcopy(filename), ios, (typemax(Int64), 0))
  end
end

"""
    COMMENT(self::AbaqusExporter, Text::AbstractString)

Write out the `**` comment option.
"""
function COMMENT(self::AbaqusExporter, Text::AbstractString)
  println(self.ios, "**" * Text)
end

"""
    HEADING(self::AbaqusExporter, Text::AbstractString)

Write out the `*HEADING` option.
"""
function HEADING(self::AbaqusExporter, Text::AbstractString)
  println(self.ios, "*HEADING");
  println(self.ios, Text)
end

"""
    PART(self::AbaqusExporter, NAME::AbstractString)

Write out the `*PART` option.
"""
function PART(self::AbaqusExporter, NAME::AbstractString)
  println(self.ios, "*PART, NAME=" * NAME)
end

"""
    END_PART(self::AbaqusExporter)

Write out the `*END PART` option.
"""
function END_PART(self::AbaqusExporter)
  println(self.ios, "*END PART")
end

"""
    ASSEMBLY(self::AbaqusExporter, NAME::AbstractString)

Write out the `*ASSEMBLY` option.
"""
function ASSEMBLY(self::AbaqusExporter, NAME::AbstractString)
  println(self.ios, "*ASSEMBLY, NAME=" * NAME)
end

"""
    END_ASSEMBLY(self::AbaqusExporter)

Write out the `*END ASSEMBLY` option.
"""
function END_ASSEMBLY(self::AbaqusExporter)
  println(self.ios, "*END ASSEMBLY")
end

"""
    INSTANCE(self::AbaqusExporter, NAME::AbstractString, PART::AbstractString)

Write out the `*INSTANCE` option.
"""
function INSTANCE(self::AbaqusExporter, NAME::AbstractString, PART::AbstractString)
  println(self.ios, "*INSTANCE, NAME=" * NAME * ", PART=" * PART)
end

"""
    END_INSTANCE(self::AbaqusExporter)

Write out the `*END INSTANCE` option.
"""
function END_INSTANCE(self::AbaqusExporter)
  println(self.ios, "*END INSTANCE")
end

"""
    NODE(self::AbaqusExporter, xyz::AbstractArray{T, 2}) where {T}

Write out the `*NODE` option.

`xyz`=array of node coordinates
"""
function NODE(self::AbaqusExporter, xyz::AbstractArray{T, 2}) where {T}
  println(self.ios, "*NODE")
  for j = 1:size(xyz,1)
    print(self.ios, "$j,$(xyz[j, 1])")
    for ixxxx = 2:1:size(xyz,2)
      print(self.ios, ",$(xyz[j, ixxxx])")
    end
    print(self.ios, "\n")
  end
end

"""
    ELEMENT(self::AbaqusExporter, TYPE::AbstractString, ELSET::AbstractString,
      start::Integer, conn::AbstractArray{T, 2}) where {T<:Integer}

Write out the `*ELEMENT` option.

`TYPE`= element type code,
`ELSET`= element set to which the elements belong,
`start`= start the element numbering at this integer,
`conn`= connectivity array for the elements, one row per element
"""
function ELEMENT(self::AbaqusExporter, TYPE::AbstractString, ELSET::AbstractString,
    start::Integer, conn::AbstractArray{T, 2}) where {T<:Integer}
    # Check that start is valid
    @assert start > 0 "The  starting element number must be > 0"
    # Check that the current element number range is disjoint from the
    # range of the elements on input: they must not overlap
    @assert (self.element_range[2] < start) ||
            (self.element_range[1] > start+size(conn, 1)-1) "Elements  must be given unique numbers"
    # Update the element range
    self.element_range = (min(self.element_range[1], start),
                            max(self.element_range[2], start+size(conn, 1)-1))
    println(self.ios, "*ELEMENT, TYPE =" * TYPE * ", ELSET=" * ELSET)
    for j=1:size(conn,1)
        print(self.ios, "$(j+start-1),")
        for ixxxx = 1:size(conn,2)-1
        if ixxxx>15
            print(self.ios, "\n")
        end
        print(self.ios, "$(conn[j,ixxxx]),")
        end
        print(self.ios, "$(conn[j,size(conn,2)])\n")
    end
end

function ELEMENT(self::AbaqusExporter, TYPE::AbstractString, ELSET::AbstractString,
    conn::AbstractArray{T, 2}) where {T<:Integer}
    start = self.element_range[2] + 1
    ELEMENT(self, TYPE, ELSET, start, conn)
end

function ELEMENT(self::AbaqusExporter, TYPE::AbstractString, ELSET::AbstractString,
    start::Integer, conn::C) where {C}
    @assert typeof(conn[1])<:Tuple
    c = fill(0, length(conn), length(conn[1]))
    for i = 1:length(conn)
        c[i, :] = [conn[i]...]
    end
    ELEMENT(self, TYPE, ELSET, start, c)
end

function ELEMENT(self::AbaqusExporter, TYPE::AbstractString, ELSET::AbstractString,
    conn::C) where {C}
    start = self.element_range[2] + 1
    ELEMENT(self, TYPE, ELSET, start, c)
end

"""
    NSET_NSET(self::AbaqusExporter, NSET::AbstractString,
      n::AbstractVector{T}) where {T<:Integer}

Write out the `*NSET` option.

`NSET` = name of the set,
`n` = array of the node numbers
"""
function NSET_NSET(self::AbaqusExporter, NSET::AbstractString, n::AbstractVector{T}) where {T<:Integer}
  println(self.ios, "*NSET, NSET=" * NSET)
  for j=1:length(n)-1
    println(self.ios, "$(n[j]),")
  end
  println(self.ios, "$(n[length(n)])")
end

"""
    ELSET_ELSET(self::AbaqusExporter, ELSET::AbstractString, n::AbstractArray{T, 1}) where {T<:Integer}

Write out the `*ELSET` option.

`ELSET` = name of the set,
`n` = array of the node numbers
"""
function ELSET_ELSET(self::AbaqusExporter, ELSET::AbstractString, n::AbstractArray{T, 1}) where {T<:Integer}
  println(self.ios, "*ELSET, ELSET=" * ELSET)
  for j=1:length(n)-1
    println(self.ios, "$(n[j]),")
  end
  println(self.ios, "$(n[length(n)])")
end

"""
    ORIENTATION(self::AbaqusExporter, ORIENTATION::AbstractString,
      a::AbstractArray{T,1}, b::AbstractArray{T,1})

Write out the `*ORIENTATION` option.

Invoke at level: Part,  Part instance,  Assembly
"""
function ORIENTATION(self::AbaqusExporter, ORIENTATION::AbstractString,
  a::AbstractArray{T,1}, b::AbstractArray{T,1}) where {T<:Real}
  println(self.ios, "*ORIENTATION,NAME=" * ORIENTATION)
  println(self.ios, "$(a[1]),$(a[2]),$(a[3]),$(b[1]),$(b[2]),$(b[3])")
  println(self.ios, "1,0.0")
end

"""
    MATERIAL(self::AbaqusExporter, MATERIAL::AbstractString)

Write out the `*MATERIAL` option.
"""
function MATERIAL(self::AbaqusExporter, MATERIAL::AbstractString)
  println(self.ios, "*MATERIAL,NAME=" * MATERIAL);
end

"""
    ELASTIC(self::AbaqusExporter, E1::F, E2::F, E3::F, nu12::F, nu13::F, nu23::F,
      G12::F, G13::F, G23::F) where {F}

Write out the `*ELASTIC,TYPE=ENGINEERING CONSTANTS` option.
"""
function ELASTIC(self::AbaqusExporter, E1::F, E2::F, E3::F, nu12::F, nu13::F, nu23::F,
  G12::F, G13::F, G23::F) where {F}
  println(self.ios, "*ELASTIC,TYPE=ENGINEERING CONSTANTS");
  println(self.ios, "$E1,$E2,$E3,$nu12,$nu13,$nu23,$G12,$G13,");
  println(self.ios, "$G23");
end

"""
    ELASTIC(self::AbaqusExporter, E::F, nu::F) where {F}

Write out the `*ELASTIC,TYPE=ISOTROPIC` option.
"""
function ELASTIC(self::AbaqusExporter, E::F, nu::F) where {F}
  println(self.ios, "*ELASTIC,TYPE=ISOTROPIC ");
  println(self.ios, "$E,$nu");
end

"""
    EXPANSION(self::AbaqusExporter, CTE::F) where {F}

Write out the `*EXPANSION` option.
"""
function EXPANSION(self::AbaqusExporter, CTE::F) where {F}
  println(self.ios, "*EXPANSION");
  println(self.ios, "$CTE,");
end

"""
    DENSITY(self::AbaqusExporter, rho)

Write out the `*DENSITY` option.
"""
function DENSITY(self::AbaqusExporter, rho)
  println(self.ios, "*DENSITY ");
  println(self.ios, "$rho");
end

"""
    SECTION_CONTROLS(self::AbaqusExporter, NAME::AbstractString,
      OPTIONAL::AbstractString)

Write out the `*SECTION CONTROLS` option.

`OPTIONAL` = string, for instance
                     HOURGLASS=ENHANCED
"""
function SECTION_CONTROLS(self::AbaqusExporter, NAME::AbstractString,
  OPTIONAL::AbstractString)
  println(self.ios, "*SECTION CONTROLS, NAME=" * NAME * "," * OPTIONAL);
end

"""
    SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
      ORIENTATION::AbstractString, ELSET::AbstractString,
      CONTROLS::AbstractString)

Write out the `*SOLID SECTION` option.

Level: Part,  Part instance

"""
function SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
  ORIENTATION::AbstractString, ELSET::AbstractString,
  CONTROLS::AbstractString)
  println(self.ios, "*SOLID SECTION,MATERIAL=" * MATERIAL *
      ",ORIENTATION =" * ORIENTATION * ",ELSET=" * ELSET *
      ", CONTROLS  =" * CONTROLS);
end

"""
    SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
      ORIENTATION::AbstractString, ELSET::AbstractString)

Write out the `*SOLID SECTION` option.

Level: Part,  Part instance

"""
function SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
  ORIENTATION::AbstractString, ELSET::AbstractString)
  println(self.ios, "*SOLID SECTION,MATERIAL=" * MATERIAL *
      ",ORIENTATION =" * ORIENTATION * ",ELSET=" * ELSET);
end

"""
    SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
      ORIENTATION::AbstractString, ELSET::AbstractString)

Write out the `*SOLID SECTION` option.

Level: Part,  Part instance

"""
function SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
  ORIENTATION::AbstractString, ELSET::AbstractString, thickness::F) where {F}
  SOLID_SECTION(self, MATERIAL, ORIENTATION, ELSET)
  println(self.ios, "$(thickness),");
end

# """
#     HOURGLASS(self::AbaqusExporter, KIND::AbstractString, VALUE::F) where {F}
#
# Write out the `*HOURGLASS` option.
#
# ```
# Example:
# *SOLID SECTION,ELSET=SOLID3,MATERIAL=MAT,CONTROL=Ac
# *HOURGLASS STIFFNESS
# 5.E8
# ```
# """
# function HOURGLASS(self::AbaqusExporter, KIND::AbstractString, VALUE::F) where {F}
#   println(self.ios, "*HOURGLASS " * KIND);
#   println(self.ios, "$( VALUE )");
# end

"""
    SURFACE_SECTION(self::AbaqusExporter, ELSET::AbstractString)

Write out the `*SURFACE SECTION` option.
"""
function SURFACE_SECTION(self::AbaqusExporter, ELSET::AbstractString)
  println(self.ios, "*SURFACE SECTION, ELSET=" * ELSET);
end

"""
    STEP_PERTURBATION_STATIC(self::AbaqusExporter)

Write out the `*STEP,PERTURBATION` option for linear static analysis.
"""
function STEP_PERTURBATION_STATIC(self::AbaqusExporter)
  println(self.ios, "*STEP,PERTURBATION");
  println(self.ios, "*STATIC");
end

"""
    STEP_PERTURBATION_BUCKLE(self::AbaqusExporter, neigv::Integer)

Write out the `*STEP,PERTURBATION` option for linear buckling analysis.
"""
function STEP_PERTURBATION_BUCKLE(self::AbaqusExporter, neigv::Integer)
  println(self.ios, "*STEP, name=Buckling, nlgeom=NO, perturbation");
  println(self.ios, "*BUCKLE");
  println(self.ios, "$(neigv), , , , ");
end

"""
    STEP_FREQUENCY(self::AbaqusExporter, Nmodes::Integer)

Write out the `*STEP,FREQUENCY` option.
"""
function STEP_FREQUENCY(self::AbaqusExporter, Nmodes::Integer)
  println(self.ios, "*STEP");
  println(self.ios, "*FREQUENCY, EIGENSOLVER=LANCZOS");
  println(self.ios, "$(Nmodes), , ,-1.E6 \n");
end

"""
    BOUNDARY(self::AbaqusExporter, NSET::AbstractString,
        is_fixed::AbstractArray{B,2},  fixed_value::AbstractArray{F,2}) where {B, F}

Write out the `*BOUNDARY` option.

`is_fixed`= array of Boolean flags (true means fixed, or prescribed),  one row per node,
`fixed_value`=array of displacements to which the corresponding displacement components is fixed
"""
function BOUNDARY(self::AbaqusExporter, NSET::AbstractString,
    is_fixed::AbstractArray{B,2},  fixed_value::AbstractArray{F,2}) where {B, F}
  println(self.ios, "*BOUNDARY");
  for j=1:size(is_fixed,1)
    for k=1:size(is_fixed,2)
      #<node number>, <first dof>, <last dof>, <magnitude of displacement>
      if is_fixed[j,k]
        println(self.ios, NSET * "$j,$k,$k,$(fixed_value[j,k])");
      end
    end
  end
end

"""
    BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer)

Write out the `*BOUNDARY` option to fix displacements at zero for a node set.

Invoke at Level: Model,  Step

`NSET`= node set,
`dof`=Degree of freedom, 1, 2, 3
"""
function BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer)
  println(self.ios, "*BOUNDARY");
  println(self.ios, NSET * ",$dof");
end

"""
    BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer,
      value::F) where {F}

Write out the `*BOUNDARY` option to fix displacements at nonzero value for a
node set.

`NSET`= node set,
`dof`=Degree of freedom, 1, 2, 3
"""
function BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer, value::F) where {F}
  println(self.ios, "*BOUNDARY,TYPE=DISPLACEMENT");
  println(self.ios, NSET * ",$dof,$dof,$value");
end

"""
    DLOAD(self::AbaqusExporter, ELSET::AbstractString,
      traction::AbstractVector{F}) where {F}

Write out the `*DLOAD` option.
"""
function DLOAD(self::AbaqusExporter, ELSET::AbstractString,
  traction::AbstractVector{F}) where {F}
  println(self.ios, "*DLOAD, follower=NO");
  nt = norm(traction)
  print(self.ios, ELSET * ",TRVEC," * "$(nt),")
  println(self.ios, "$(traction[1]/nt),$(traction[2]/nt),$(traction[3]/nt)");
end

"""
    CLOAD(self::AbaqusExporter, NSET::AbstractString, dof::Integer,
      magnitude::F) where {F}

Write out the `*CLOAD` option.

NSET=Node set
dof= 1, 2, 3,
magnitude= signed multiplier
"""
function CLOAD(self::AbaqusExporter, NSET::AbstractString, dof::Integer,
  magnitude::F) where {F}
  println(self.ios, "*CLOAD");
  println(self.ios, NSET * ",$dof,$magnitude");
end

"""
    CLOAD(self::AbaqusExporter, nodenumber::Integer, dof::Integer,
      magnitude::F) where {F}

Write out the `*CLOAD` option.

nodenumber=Number of node
dof= 1, 2, 3,
magnitude= signed multiplier
"""
function CLOAD(self::AbaqusExporter, nodenumber::Integer, dof::Integer,
  magnitude::F) where {F}
  println(self.ios, "*CLOAD");
  println(self.ios, "$(nodenumber),$dof,$magnitude");
end

"""
    TEMPERATURE(self::AbaqusExporter, nlist::AbstractArray{I, 1},
      tlist::AbstractArray{F, 1}) where {I, F}

Write out the `*TEMPERATURE` option.
"""
function TEMPERATURE(self::AbaqusExporter, Classifier::AbstractString,
  nlist::AbstractArray{I, 1},
  tlist::AbstractArray{F, 1}) where {I, F}
  @assert length(nlist) == length(tlist)
  println(self.ios, "*TEMPERATURE");
  for ixxxx = 1:length(nlist)
    println(self.ios, Classifier * "$(nlist[ixxxx]),$(tlist[ixxxx])");
  end
end

"""
    NODE_PRINT(self::AbaqusExporter, NSET::AbstractString)

Write out the `*NODE PRINT` option.
"""
function NODE_PRINT(self::AbaqusExporter, NSET::AbstractString; KEYS = "U")
  println(self.ios, "*NODE PRINT, NSET=" * NSET);
  println(self.ios, KEYS);
end

"""
    EL_PRINT(self::AbaqusExporter, ELSET::AbstractString, KEYS::AbstractString)

Write out the `*EL PRINT` option.
"""
function EL_PRINT(self::AbaqusExporter, ELSET::AbstractString, KEYS::AbstractString)
  println(self.ios, "*EL PRINT, ELSET=" * ELSET * ", POSITION=INTEGRATION POINTS, SUMMARY= YES");
  println(self.ios, "");
  println(self.ios, KEYS);
end

"""
    ENERGY_PRINT(self::AbaqusExporter)

Write out the `*ENERGY PRINT` option.
"""
function ENERGY_PRINT(self::AbaqusExporter)
  println(self.ios, "*ENERGY PRINT");
end

"""
    END_STEP(self::AbaqusExporter)

Write out the `*END STEP` option.
"""
function END_STEP(self::AbaqusExporter)
  println(self.ios, "*END STEP");
end

"""
    close(self::AbaqusExporter)

Close  the stream opened by the exporter.
"""
function close(self::AbaqusExporter)
  close(self.ios)
end


################################################################################
# CSV export
################################################################################

"""
    savecsv(name::String; kwargs...)

Save arrays as a CSV file.
"""
function savecsv(name::String; kwargs...)
    colnames = Symbol[]
    columns = Any[]
    for apair in pairs(kwargs)
        sy, val = apair
        push!(colnames, sy)
        push!(columns, val)
    end
    ncol = length(colnames)
    nrow = length(columns[1])
    for j = 1:ncol
        @assert length(columns[j]) == nrow "Columns must have the same number of rows"
    end
    if !occursin(r"^.*\.csv$", name) && !occursin(r"^.*\.CSV$", name)
        name = name * ".csv"
    end
    open(name, "w") do fid
        for (j, c) in enumerate(colnames)
            print(fid, "\"", c, "\"")
            j < ncol && print(fid, ",")
        end
        print(fid, "\n")
        for r = 1:nrow
            for j = 1:ncol
                print(fid, columns[j][r])
                j < ncol && print(fid, ",")
            end
            print(fid, "\n")
        end
        print(fid, "\n")
    end
    return true
end

################################################################################
# NASTRAN export
################################################################################

"""
NASTRANExporter

Export mesh to Abaqus.
"""
mutable struct NASTRANExporter
    filename::AbstractString
    ios::IO
    element_range::Tuple{Int64, Int64}
    function NASTRANExporter(filename::AbstractString)
        if match(r".*\.nas$", filename) == nothing
            filename = filename * ".nas"
        end
        ios = open(filename,"w+")
        return new(deepcopy(filename), ios, (typemax(Int64), 0))
    end
end

NASTRANtypemap = Dict{DataType, AbstractString}(FESetT4=>"CTETRA", FESetT10=>"CTETRA")

"""
    CEND(self::NASTRANExporter)

Terminate the Executive Control section.
"""
function CEND(self::NASTRANExporter)
    println(self.ios, "CEND")
end

"""
    BEGIN_BULK(self::NASTRANExporter)

Terminate the Case Control section by starting the bulk section.
"""
function BEGIN_BULK(self::NASTRANExporter)
    println(self.ios, "BEGIN BULK")
end

"""
    ENDDATA(self::NASTRANExporter)

Terminate the bulk section.
"""
function ENDDATA(self::NASTRANExporter)
    println(self.ios, "ENDDATA")
end

"""
    GRID(self::NASTRANExporter, n::Int, xyz::Vector{FFlt})

Write a grid-point statement.
"""
function GRID(self::NASTRANExporter, n::Int, xyz::Vector{FFlt})
    @printf self.ios "GRID,%d,,%g,%g,%g\n" n xyz[1] xyz[2] xyz[3]
end

"""
    PSOLID(self::NASTRANExporter, pid::Int, mid::Int)

Write solid-property statement.
"""
function PSOLID(self::NASTRANExporter, pid::Int, mid::Int)
    @printf self.ios "PSOLID,%d,%d\n" pid mid
end

"""
    MAT1(self::NASTRANExporter, mid::Int, E::FFlt, G::FFlt, nu::FFlt, rho::FFlt, A::FFlt, TREF::FFlt, GE::FFlt)

Write a statement for an isotropic elastic material.
"""
function MAT1(self::NASTRANExporter, mid::Int, E::FFlt, G::FFlt, nu::FFlt, rho::FFlt, A::FFlt, TREF::FFlt, GE::FFlt)
    @printf self.ios "MAT1,%d,%g,%g,%g,%g,%g,%g,%g\n" mid E G nu rho A TREF GE
end

"""
    MAT1(self::NASTRANExporter, mid::Int, E::FFlt, nu::FFlt, rho::FFlt = 0.0, A::FFlt = 0.0, TREF::FFlt = 0.0)

Write a statement for an isotropic elastic material.
"""
function MAT1(self::NASTRANExporter, mid::Int, E::FFlt, nu::FFlt, rho::FFlt = 0.0, A::FFlt = 0.0, TREF::FFlt = 0.0, GE::FFlt = 0.0)
    @printf self.ios "MAT1,%d,%g,,%g,%g,%g,%g,%g\n" mid E nu rho A TREF GE
end

function CTETRA(self::NASTRANExporter, eid::Int, pid::Int, conn::Vector{Int})
    nc = length(conn)
    @printf self.ios "CTETRA,%d,%d" eid pid
    for k = 1:nc
        @printf self.ios ",%d" conn[k]
        if k == 6
            @printf self.ios "\n"
        end
    end
    @printf self.ios "\n"
end

"""
    close(self::NASTRANExporter)

Close  the stream opened by the exporter.
"""
function close(self::NASTRANExporter)
    close(self.ios)
end


################################################################################
# STL export
################################################################################

"""
STLExporter

Export surface mesh as STL file.
"""
mutable struct STLExporter
    filename::AbstractString
    ios::IO
    element_range::Tuple{Int64, Int64}
    function STLExporter(filename::AbstractString)
        if match(r".*\.stl$", filename) == nothing
            filename = filename * ".stl"
        end
        ios = open(filename,"w+")
        return new(deepcopy(filename), ios, (typemax(Int64), 0))
    end
end

function solid(self::STLExporter, name::AbstractString = "thesolid")
    @printf self.ios "solid %s\n" name
end

function facet(self::STLExporter, v1::Vector{FFlt}, v2::Vector{FFlt}, v3::Vector{FFlt})
    V = v2 - v1
    W = v3 - v1
    normal = cross(V, W)
    normal = normal / norm(normal)
    @printf self.ios "facet normal %e %e %e\n" normal[1] normal[2] normal[3]
    @printf self.ios "    outer loop\n"
    @printf self.ios "        vertex %e %e %e\n" v1[1] v1[2] v1[3]
    @printf self.ios "        vertex %e %e %e\n" v2[1] v2[2] v2[3]
    @printf self.ios "        vertex %e %e %e\n" v3[1] v3[2] v3[3]
    @printf self.ios "    endloop\n"
    @printf self.ios "endfacet\n"
end

function endsolid(self::STLExporter, name::AbstractString = "thesolid")
    @printf self.ios "endsolid %s\n" name
end

"""
    close(self::STLExporter)

Close  the stream opened by the exporter.
"""
function close(self::STLExporter)
    close(self.ios)
end


################################################################################
# H2Lib triangular-surface export
################################################################################

import Base.BitSet

mutable struct _HypFace
    n::Int # serial number of the hyperface
    hf::Array{Int} # numbers of vertices defining the hyperface
    a::Int # anchor node number
    o::BitSet # numbers of the other nodes on the hyperface
end

mutable struct _HypFaceContainer
    nv::Int # total number of vertices per hyperface
    d::Dict{Int, Array{_HypFace}} # dictionary
    n::Int # total number of hyperfaces in the container
end

makehypfacedict(nv) = _HypFaceContainer(nv, Dict{Int, Array{_HypFace}}(), 0);

function hypfacen(container, hypf)
    h = sort([i for i in hypf])
    anchor = h[1]; other = BitSet(h[2:end]);
    ca = get(container.d, anchor, _HypFace[]);
    for k = 1:length(ca)
        if (ca[k].a == anchor) && (ca[k].o == other)
            return ca[k].n # do have this hyper face: here is its serial number
        end
    end
    return 0 # do not have this hyper face
end

function inserthypface!(container, hypf)
    if hypfacen(container, hypf) == 0
        @assert length(hypf) == container.nv "Hyperface does not have the right number of vertices"
        h = sort([i for i in hypf])
        anchor = h[1]; other = BitSet(h[2:end]);
        ca = get(container.d, anchor, _HypFace[]);
        container.n = container.n + 1
        push!(ca, _HypFace(container.n, hypf, anchor, other));
        container.d[anchor] = ca
    end
    return container
end

function hypfacedicttoarray(container)
    a = fill(0, container.n, container.nv)
    for ca in values(container.d)
        for i = 1:length(ca)
            a[ca[i].n, :] = ca[i].hf
        end
    end
    return a
end

function h2libexporttri(theFile::String, Connectivity, Points)
    @assert size(Connectivity, 2) == 3 "Only triangles accepted"
    adjust = 1 # the library is C-based and expects indexes to be zero-based

    ed = makehypfacedict(2)
    for i = 1:size(Connectivity, 1)
        inserthypface!(ed, Connectivity[i, [1, 2]])
        inserthypface!(ed, Connectivity[i, [2, 3]])
        inserthypface!(ed, Connectivity[i, [3, 1]])
    end
    ea = hypfacedicttoarray(ed)

    fid=open(theFile,"w");
    if (fid == -1)
        error(["Could not open " * theFile])
        return nothing
    end

    print(fid, size(Points, 1), " ", size(ea, 1), " ", size(Connectivity, 1), "\n");

    for i= 1:size(Points, 1)
        for j= 1:size(Points,2)-1
            print(fid,Points[i,j]," ");
        end
        print(fid,Points[i,end],"\n");
    end

    for i= 1:size(ea, 1)
        for j= 1:size(ea,2)-1
            print(fid, ea[i,j]-adjust, " ");
        end
        print(fid, ea[i,end]-adjust, "\n");
    end

    for i= 1:size(Connectivity, 1)
        s = [hypfacen(ed, Connectivity[i, [2, 3]]) hypfacen(ed, Connectivity[i, [3, 1]]) hypfacen(ed, Connectivity[i, [1, 2]])]
        for j= 1:size(Connectivity,2)
            print(fid, Connectivity[i,j]-adjust, " ");
        end
        for j= 1:length(s)-1
            print(fid, s[j]-adjust, " ");
        end
        print(fid, s[end]-adjust, "\n");
    end

    fid=close(fid);
    return true
end

end
