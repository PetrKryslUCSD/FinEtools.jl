"""
    FieldModule

Module for abstract fields.
"""
module FieldModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import Base.copyto!

"""
    Abstract field.

Expected  attributes:
  + values::FMat{T} # Array of degree of freedom parameters,  indexed by node
  + dofnums::FIntMat # Array of degree of freedom numbers, indexed by node
  + is_fixed::Matrix{Bool} # Array of Boolean flags, indexed by node
  + fixed_values::FMat{T} # Array of fixed values, indexed by node
  + nfreedofs::FInt # Total number of free degrees of freedom
See the macro `add_Field_fields()`.
"""
abstract type Field end


macro add_Field_fields()
    return esc(:(
    values::Array{T, 2};
    dofnums::Array{FInt, 2};
    is_fixed::Matrix{Bool};
    fixed_values::Array{T, 2};;
    nfreedofs::FInt;
    )
    )
end

"""
    ndofs{F<:Field}(self::T)

Dimension of the degree of freedom parameters (i. e.  how many degrees of
freedom per entity).
"""
ndofs(self::Field)  = size(self.values, 2)

"""
    nents{F<:Field}(self::F)

Number of nodes associated with the field.
"""
nents(self::Field)  = size(self.values, 1)

"""
    copyto!(DEST::F,  SRC::F) where {F<:Field}

Copy data from one field to another.
"""
function copyto!(DEST::F,  SRC::F) where {F<:Field}
    copyto!(DEST.values, SRC.values)
    copyto!(DEST.dofnums, SRC.dofnums)
    copyto!(DEST.is_fixed, SRC.is_fixed)
    copyto!(DEST.fixed_values, SRC.fixed_values)
    DEST.nfreedofs = SRC.nfreedofs
    return  DEST
end

"""
    wipe!(self::Field)

Wipe all the data from the field.
"""
function wipe!(self::Field)
    Zer = zero(eltype(self.fixed_values[1]))
    self.nfreedofs = 0
    fill!(self.dofnums, 0)
    fill!(self.is_fixed, false)
    fill!(self.values, Zer)
    fill!(self.fixed_values, Zer)
    return self
end

"""
    gathersysvec{F<:Field}(self::F)

Gather values from the field for the whole system vector.
"""
function gathersysvec(self::F) where {F<:Field}
    nents,dim = size(self.values)
    vec = zeros(typeof(self.values[1,1]),self.nfreedofs,1)
    for i = 1:nents
        for j = 1:dim
            en = self.dofnums[i,j]
            if (en > 0) && (en <= self.nfreedofs)
                vec[en] = self.values[i,j]
            end
        end
    end
    return vec
end

"""
    gathersysvec!(self::F, vec::FVec{T}) where {F<:Field, T}

Gather values from the field for the whole system vector.
"""
function gathersysvec!(self::F, vec::FVec{T}) where {F<:Field, T}
    nents,dim = size(self.values)
    @assert length(vec) == self.nfreedofs
    for i = 1:nents
        for j = 1:dim
            en = self.dofnums[i,j]
            if (en > 0) && (en <= self.nfreedofs)
                vec[en] = self.values[i,j]
            end
        end
    end
    return vec
end

"""
    gathervalues_asvec!(self::Field, dest::AbstractArray{T, 1},
        conn::CC) where {CC, T}

Gather values from the field into a vector.

The order is: for each node  in the connectivity, copy into the buffer all the
degrees of freedom,  then the next node and so on.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gathervalues_asvec!(self::Field, dest::AbstractArray{T, 1}, conn::CC) where {CC, T}
    en::FInt = 1;
    for i = 1:length(conn)
        for j = 1:size(self.values,2)
            dest[en] = self.values[conn[i],j];
            en = en + 1;
        end
    end
    return dest
end

"""
    gathervalues_asmat!(self::Field, dest::AbstractArray{T, 2},
        conn::CC) where {CC, T}

Gather values from the field into a two-dimensional array.

The order is: for each node  in the connectivity, copy into the corresponding
row of the buffer all the degrees of freedom,  then the next node into the next
row and so on.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gathervalues_asmat!(self::Field, dest::AbstractArray{T, 2},    conn::CC) where {CC, T}
    for i = 1:length(conn)
        for j = 1:size(self.values,2)
            dest[i, j] = self.values[conn[i], j];
        end
    end
    return dest
end

"""
    gatherfixedvalues_asvec!(self::Field, dest::AbstractArray{T, 1},
        conn::CC) where {CC, T}

Gather FIXED values from the field into a vector.

The order is: for each node  in the connectivity, copy into the buffer all the
fixed degrees of freedom,  then the next node and so on. If a degree of freedom
is NOT fixed, the corresponding entry is  set to zero.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gatherfixedvalues_asvec!(self::Field, dest::AbstractArray{T, 1},    conn::CC) where {CC, T}
    en::FInt = 1;
    for i = 1:length(conn)
        for j = 1:size(self.fixed_values,2)
            if self.is_fixed[conn[i],j] # free degree of freedom
                dest[en] = self.fixed_values[conn[i], j];
            else
                dest[en] = 0.0
            end
            en = en + 1;
        end
    end
    return dest
end

"""
    gatherfixedvalues_asmat!(self::Field, dest::AbstractArray{T, 2},
        conn::CC) where {CC, T}

Gather FIXED values from the field into a two-dimensional array.

The order is: for each node  in the connectivity, copy into the corresponding
row of the buffer all the degrees of freedom,  then the next node into the next
row and so on.  If a degree of freedom
is NOT fixed, the corresponding entry is  set to zero.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gatherfixedvalues_asmat!(self::Field, dest::AbstractArray{T, 2},    conn::CC) where {CC, T}
    for i = 1:length(conn)
        for j = 1:size(self.fixed_values,2)
            if self.is_fixed[conn[i],j] # fixed degree of freedom
                dest[i, j] = self.fixed_values[conn[i], j];
            else
                dest[i, j] = 0.0
            end
        end
    end
    return dest
end

function anyfixedvaluenz(self::Field, conn::CC) where {CC}
    for i = 1:length(conn)
        for j = 1:size(self.fixed_values,2)
            if self.is_fixed[conn[i],j] # free degree of freedom
                if  abs(self.fixed_values[conn[i], j]) > 0.0
                    return true
                end 
            end
        end
    end
    return false
end

"""
    gatherdofnums!(self::Field, dest::A, conn::CC) where {A, CC}

Gather dofnums from the field.
"""
function gatherdofnums!(self::Field, dest::A, conn::CC) where {A, CC}
    en::FInt = 1;
    for i = 1:length(conn)
        for j = 1:size(self.dofnums,2)
            dest[en] = self.dofnums[conn[i],j];
            en = en+1;
        end
    end
    return dest
end

"""
    numberdofs!(self::Field)

Number the degrees of freedom.

The free components in the field are numbered consecutively. No effort is
made to optimize the numbering in any way. If you'd like to optimize the
numbering of the degrees of freedom, use the above form that sets the
permutation of the degrees of freedom, or the permutation of the nodes.
"""
function numberdofs!(self::Field)
    fixed_dofnum::FInt = 0
    nents,dim = size(self.values)
    self.nfreedofs::FInt =0
    for i=1:nents
        for j=1:dim
            if !self.is_fixed[i,j] # free degree of freedom
                self.nfreedofs = self.nfreedofs + 1
                self.dofnums[i,j] = self.nfreedofs
            else # fixed degree of freedom: no equation
                self.dofnums[i,j] = fixed_dofnum
            end
        end
    end
    return  self
end

"""
    setebc!(self::Field, fenids::FIntVec, is_fixed::Bool, comp::FInt,
      val::FVec{T}) where {T<:Number}

Set the EBCs (essential boundary conditions).

`fenids`         - array of N node identifiers
`is_fixed` = scaler Boolean: are the degrees of freedom being fixed (true)
             or released (false),
`comp` = integer, which  degree of freedom (component),
`val` = array of N values of type T

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::Field, fenids::FIntVec, is_fixed::Bool, comp::FInt,    val::FVec{T}) where {T<:Number}
    @assert comp <= size(self.values,2) "Requested  nonexistent  degree of freedom"
    @assert maximum(fenids) <= size(self.values,1) "Requested nonexistent node"
    @assert size(fenids) == size(val) "Arrays of mismatched sizes"
    for  j = 1:length(fenids)
        self.is_fixed[fenids[j],comp] = is_fixed;
        if self.is_fixed[fenids[j],comp]
            self.fixed_values[fenids[j],comp] = val[j];
        else
            self.fixed_values[fenids[j],comp] = zero(T)
        end
    end
    self.nfreedofs = 0
    fill!(self.dofnums, 0)
    return  self
end

"""
    setebc!(self::Field, fenids::FIntVec, is_fixed::Bool, comp::FInt,
      val::T) where {T<:Number}

Set the EBCs (essential boundary conditions).

`fenids`         - array of N node identifiers
`is_fixed` = scaler Boolean: are the degrees of freedom being fixed (true)
             or released (false),
`comp` = integer, which  degree of freedom (component),
`val` = scalar of type T

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::Field, fenids::FIntVec, is_fixed::Bool, comp::FInt,    val::T) where {T<:Number}
    @assert (comp >= 1 && comp <= size(self.values,2)) "Requested  nonexistent  degree of freedom"
    @assert maximum(fenids) <= size(self.values,1) "Requested nonexistent node"
    @assert minimum(fenids) >= 1 "Requested nonexistent node"
    for  j = 1:length(fenids)
        self.is_fixed[fenids[j],comp] = is_fixed;
        if self.is_fixed[fenids[j],comp]
            self.fixed_values[fenids[j],comp] = val;
        else
            self.fixed_values[fenids[j],comp] = zero(T)
        end
    end
    self.nfreedofs = 0
    fill!(self.dofnums, 0)
    return  self
end

"""
    setebc!(self::Field, fenids::FIntVec, comp::FInt,
      val::FVec{T}) where {T<:Number}

Set the EBCs (essential boundary conditions).

`fenids`         - array of N node identifiers
`comp` = integer, which  degree of freedom (component),
`val` = array of N values of type T

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::Field, fenids::FIntVec, comp::FInt, val::FVec{T}) where {T<:Number}
    return setebc!(self, fenids, true, comp, val)
end


"""
    setebc!(self::Field, fenids::FIntVec, comp::FInt;
      val::T=0.0) where {T<:Number}

Set the EBCs (essential boundary conditions).

`fenids`         - array of N node identifiers
`comp` = integer, which  degree of freedom (component),
`val` = scalar of type T

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::Field, fenids::FIntVec, comp::FInt; val::T=0.0) where {T<:Number}
    return setebc!(self, fenids, true, comp, val)
end

"""
    setebc!(self::Field, fenids::FIntVec, is_fixed::Bool, comp::FInt;
      val::T=0.0) where {T<:Number}

Set the EBCs (essential boundary conditions).

`fenids`         - array of N node identifiers
`comp` = integer, which  degree of freedom (component),
`val` = scalar of type T

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::Field, fenids::FIntVec, is_fixed::Bool, comp::FInt;    val::T=0.0) where {T<:Number}
    j = comp
    @assert (j >= 1) && (j <= ndofs(self))
    setebc!(self, fenids, is_fixed, j, val)
    return self
end

"""
    setebc!(self::Field, fenids::FIntVec)

Set the EBCs (essential boundary conditions).

Suppress all degrees of freedom at the given nodes.

`fenids`         - array of N node identifiers

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::Field, fenids::FIntVec)
    Zer = zero(eltype(self.fixed_values[1]))
    for comp = 1:size(self.values, 2)
        setebc!(self, fenids, true, comp, Zer)
    end
    return self
end

"""
    setebc!(self::Field, fenid::FInt)

Set the EBCs (essential boundary conditions).

Suppress all degrees of freedom at the given node.

`fenid`         - One integer as a node identifier

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::Field, fenid::FInt)
    Zer = zero(eltype(self.fixed_values[1]))
    for comp = 1:size(self.values, 2)
        setebc!(self, [fenid], true, comp, Zer)
    end
    return self
end

"""
    setebc!(self::Field)

Set the EBCs (essential boundary conditions).

All essential boundary conditions are CLEARED.

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::Field)
    self.nfreedofs = 0
    fill!(self.dofnums, 0)
    fill!(self.is_fixed, false)
    fill!(self.fixed_values, zero(eltype(self.fixed_values[1])))
    return  self
end

"""
    applyebc!(self::Field)

Apply EBCs (essential boundary conditions).
"""
function applyebc!(self::Field)
    nents,dim = size(self.values);
    for i = 1:nents
        for j = 1:dim
            if self.is_fixed[i,j]
                self.values[i,j] = self.fixed_values[i,j];
            end
        end
    end
    return  self
end

"""
    scattersysvec!(self::Field, vec::FVec{T}) where {T<:Number}

Scatter values to the field from a system vector.
"""
function scattersysvec!(self::Field, vec::FVec{T}) where {T<:Number}
    nents,dim = size(self.values);
    for i = 1:nents
        for j = 1:dim
            dn = self.dofnums[i,j];
            if (dn > 0) && (dn <= self.nfreedofs)
                self.values[i,j] = vec[dn];
            end
        end
    end
    return  self
end


end
