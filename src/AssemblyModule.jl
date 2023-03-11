"""
    AssemblyModule

Module for assemblers  of system matrices and vectors.
"""
module AssemblyModule

__precompile__(true)

using ..FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using SparseArrays: sparse, spzeros, SparseMatrixCSC
using LinearAlgebra: diag

"""
    AbstractSysmatAssembler

Abstract type of system-matrix assembler.
"""
abstract type AbstractSysmatAssembler end

"""
    SysmatAssemblerSparse{IT, MBT, IBT} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices.

!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix`.
"""
mutable struct SysmatAssemblerSparse{IT, MBT, IBT} <: AbstractSysmatAssembler
    buffer_length::IT
    matbuffer::MBT
    rowbuffer::IBT
    colbuffer::IBT
    buffer_pointer::IT
    ndofs_row::IT
    ndofs_col::IT
    nomatrixresult::Bool
end

"""
    SysmatAssemblerSparse(z= zero(FFlt), nomatrixresult = false)

Construct blank system matrix assembler. 

The matrix entries are of type `T`. The assembler either produces a sparse
matrix (when `nomatrixresult = true`), or does not (when `nomatrixresult =
false`). When the assembler does not produce the sparse matrix when
`makematrix!` is called, it still can be constructed from the buffers stored in
the assembler.


# Example

This is how a sparse matrix is assembled from two rectangular dense matrices.
```
    a = SysmatAssemblerSparse(0.0)                                                        
    startassembly!(a, 5, 5, 3, 7, 7)    
    m = [0.24406   0.599773    0.833404  0.0420141                                             
        0.786024  0.00206713  0.995379  0.780298                                              
        0.845816  0.198459    0.355149  0.224996]                                     
    assemble!(a, m, [1 7 5], [5 2 1 4])        
    m = [0.146618  0.53471   0.614342    0.737833                                              
         0.479719  0.41354   0.00760941  0.836455                                              
         0.254868  0.476189  0.460794    0.00919633                                            
         0.159064  0.261821  0.317078    0.77646                                               
         0.643538  0.429817  0.59788     0.958909]                                   
    assemble!(a, m, [2 3 1 7 5], [6 7 3 4])                                        
    A = makematrix!(a) 
```

When the `nomatrixresult` is set as true, no matrix is produced.
```
    a = SysmatAssemblerSparse(0.0, true)                                                        
    startassembly!(a, 5, 5, 3, 7, 7)    
    m = [0.24406   0.599773    0.833404  0.0420141                                             
        0.786024  0.00206713  0.995379  0.780298                                              
        0.845816  0.198459    0.355149  0.224996]                                     
    assemble!(a, m, [1 7 5], [5 2 1 4])        
    m = [0.146618  0.53471   0.614342    0.737833                                              
         0.479719  0.41354   0.00760941  0.836455                                              
         0.254868  0.476189  0.460794    0.00919633                                            
         0.159064  0.261821  0.317078    0.77646                                               
         0.643538  0.429817  0.59788     0.958909]                                   
    assemble!(a, m, [2 3 1 7 5], [6 7 3 4])                                        
    A = makematrix!(a) 
```
Here `A` is a sparse zero matrix. To construct the correct matrix is still 
possible, for instance like this:
```
    a.nomatrixresult = false
    A = makematrix!(a) 
```
At this point all the buffers of the assembler have been cleared, and 
`makematrix!(a) ` is no longer possible.

"""
function SysmatAssemblerSparse(z= zero(FFlt), nomatrixresult = false)
    return SysmatAssemblerSparse(0, FFlt[z], FInt[0], FInt[0], 0, 0, 0, nomatrixresult)
end

"""
    startassembly!(
        self::SysmatAssemblerSparse,
        elem_mat_nrows,
        elem_mat_ncols,
        elem_mat_nmatrices,
        ndofs_row,
        ndofs_col,
    )

Start the assembly of a global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.
- `elem_mat_nrows`= number of rows in typical element matrix,
- `elem_mat_ncols`= number of columns in a typical element matrix,
- `elem_mat_nmatrices`= number of element matrices,
- `ndofs_row`= Total number of equations in the row direction,
- `ndofs_col`= Total number of equations in the column direction.

If the `buffer_pointer` field of the assembler is 0, which is the case after
that assembler was created, the buffers are resized appropriately given the
dimensions on input. Otherwise, the buffers are left completely untouched.

# Returns
- `self`: the modified assembler.


!!! note

    The buffers are initially not filled with anything meaningful.
    After the assembly, only the `(self.buffer_pointer - 1)` entries
    are meaningful numbers. Beware!
"""
function startassembly!(
    self::SysmatAssemblerSparse,
    elem_mat_nrows,
    elem_mat_ncols,
    elem_mat_nmatrices,
    ndofs_row,
    ndofs_col,
)
    # Only reset the buffers if the pointer is 0.10
    if self.buffer_pointer < 1
        self.buffer_length = elem_mat_nmatrices * elem_mat_nrows * elem_mat_ncols
        resize!(self.rowbuffer, self.buffer_length)
        resize!(self.colbuffer, self.buffer_length)
        resize!(self.matbuffer, self.buffer_length)
        self.rowbuffer .= 1
        self.colbuffer .= 1
        self.matbuffer .= 0
        self.buffer_pointer = 1
        self.ndofs_row = ndofs_row
        self.ndofs_col = ndofs_col
    end
    return self
end

"""
    assemble!(
        self::SysmatAssemblerSparse,
        mat::MT,
        dofnums_row::IT,
        dofnums_col::IT,
    ) where {MT, IT}

Assemble a rectangular matrix.
"""
function assemble!(
    self::SysmatAssemblerSparse,
    mat::MT,
    dofnums_row::IT,
    dofnums_col::IT,
) where {MT, IT}
    # Assembly of a rectangular matrix.
    # The method assembles a rectangular matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    p = self.buffer_pointer
    @assert p + ncolumns * nrows <= self.buffer_length + 1
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns, i in 1:nrows
        self.matbuffer[p] = mat[i, j] # serialized matrix
        self.rowbuffer[p] = dofnums_row[i]
        self.colbuffer[p] = dofnums_col[j]
        p = p + 1
    end
    self.buffer_pointer = p
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparse)

Make a sparse matrix.
"""
function makematrix!(self::SysmatAssemblerSparse)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    @assert length(self.rowbuffer) >= self.buffer_pointer - 1
    @assert length(self.colbuffer) >= self.buffer_pointer - 1
    # Here we will go through the rows and columns, and whenever the row or
    # the column refer to indexes outside of the limits of the matrix, the
    # corresponding value will be set to 0 and assembled to row and column 1.
    @inbounds for j in 1:self.buffer_pointer-1
        if (self.rowbuffer[j] > self.ndofs_row) || (self.rowbuffer[j] <= 0)
            self.rowbuffer[j] = 1
            self.matbuffer[j] = 0.0
        end
        if (self.colbuffer[j] > self.ndofs_col) || (self.colbuffer[j] <= 0)
            self.colbuffer[j] = 1
            self.matbuffer[j] = 0.0
        end
    end
    if self.nomatrixresult
        # No actual sparse matrix is returned. The entire result of the assembly
        # is preserved in the assembler buffers. 
        return spzeros(self.ndofs_row, self.ndofs_col)
    else
        # The sparse matrix is constructed and returned. The  buffers used for
        # the assembly are cleared.
        S = sparse(
            self.rowbuffer[1:self.buffer_pointer-1],
            self.colbuffer[1:self.buffer_pointer-1],
            self.matbuffer[1:self.buffer_pointer-1],
            self.ndofs_row,
            self.ndofs_col,
        )
        self = SysmatAssemblerSparse(zero(eltype(self.matbuffer)))# get rid of the buffers
        return S
    end
end

"""
    SysmatAssemblerSparseSymm{IT, MBT, IBT} <: AbstractSysmatAssembler

Assembler for a **symmetric square** matrix  assembled from symmetric square
matrices.
"""
mutable struct SysmatAssemblerSparseSymm{IT, MBT, IBT} <: AbstractSysmatAssembler
    buffer_length::IT
    matbuffer::MBT
    rowbuffer::IBT
    colbuffer::IBT
    buffer_pointer::IT
    ndofs::IT
end

"""
    SysmatAssemblerSparseSymm(zero::T=0.0) where {T<:Number}

Construct blank system matrix assembler for symmetric matrices. The matrix
entries are of type `T`.

# Example

This is how a symmetric sparse matrix is assembled from two square dense matrices.
```
	a = SysmatAssemblerSparseSymm(0.0)                                                        
	startassembly!(a, 5, 5, 3, 7, 7)    
	m = [0.24406   0.599773    0.833404  0.0420141                                             
		0.786024  0.00206713  0.995379  0.780298                                              
		0.845816  0.198459    0.355149  0.224996]                                     
	assemble!(a, m'*m, [5 2 1 4], [5 2 1 4])        
	m = [0.146618  0.53471   0.614342    0.737833                                              
		 0.479719  0.41354   0.00760941  0.836455                                              
		 0.254868  0.476189  0.460794    0.00919633                                            
		 0.159064  0.261821  0.317078    0.77646                                               
		 0.643538  0.429817  0.59788     0.958909]                                   
	assemble!(a, m'*m, [2 3 1 5], [2 3 1 5])                                        
	A = makematrix!(a) 
```
# See also
SysmatAssemblerSparse
"""
function SysmatAssemblerSparseSymm(z= zero(FFlt))
    return SysmatAssemblerSparseSymm(0, FFlt[z], FInt[0], FInt[0], 0, 0)
end

"""
    startassembly!(
        self::SysmatAssemblerSparseSymm,
        elem_mat_dim,
        ignore1,
        elem_mat_nmatrices,
        ndofs,
        ignore2,
    )

Start the assembly of a symmetric square global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.
- `elem_mat_nrows`= number of rows in typical element matrix,
- `elem_mat_ncols`= number of columns in a typical element matrix,
- `elem_mat_nmatrices`= number of element matrices,
- `ndofs_row`= Total number of equations in the row direction,
- `ndofs_col`= Total number of equations in the column direction.

If the `buffer_pointer` field of the assembler is 0, which is the case after
that assembler was created, the buffers are resized appropriately given the
dimensions on input. Otherwise, the buffers are left completely untouched.

# Returns
- `self`: the modified assembler.


!!! note

    The buffers are initially not filled with anything meaningful.
    After the assembly, only the `(self.buffer_pointer - 1)` entries
    are meaningful numbers. Beware!
"""
function startassembly!(
    self::SysmatAssemblerSparseSymm,
    elem_mat_dim,
    ignore1,
    elem_mat_nmatrices,
    ndofs,
    ignore2,
)
    # Only reset the buffers if the pointer is 0.10
    if self.buffer_pointer < 1
        self.buffer_length = elem_mat_nmatrices * elem_mat_dim^2
        resize!(self.rowbuffer, self.buffer_length)
        resize!(self.colbuffer, self.buffer_length)
        resize!(self.matbuffer, self.buffer_length)
        self.rowbuffer .= 1
        self.colbuffer .= 1
        self.matbuffer .= 0
        self.buffer_pointer = 1
        self.ndofs = ndofs
    end
    return self
end

"""
    assemble!(
        self::SysmatAssemblerSparseSymm,
        mat::MT,
        dofnums::IT,
        ignore
    ) where {MT, IT}

Assemble a square symmetric matrix.

`dofnums` are the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are assumed to be the same).
"""
function assemble!(
    self::SysmatAssemblerSparseSymm,
    mat::MT,
    dofnums::IT,
    ignore
) where {MT, IT}
    # Assembly of a square symmetric matrix.
    # The method assembles the lower triangle of the square symmetric matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums)
    ncolumns = nrows
    p = self.buffer_pointer
    @assert p + ncolumns * nrows <= self.buffer_length + 1
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns, i in j:nrows
        self.matbuffer[p] = mat[i, j] # serialized matrix
        self.rowbuffer[p] = dofnums[i]
        self.colbuffer[p] = dofnums[j]
        p = p + 1
    end
    self.buffer_pointer = p
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparseSymm)

Make a sparse symmetric square matrix.
"""
function makematrix!(self::SysmatAssemblerSparseSymm)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    @assert length(self.rowbuffer) >= self.buffer_pointer - 1
    @assert length(self.colbuffer) >= self.buffer_pointer - 1
    # Here we will go through the rows and columns, and whenever the row or
    # the column refer to indexes outside of the limits of the matrix, the
    # corresponding value will be set to 0 and assembled to row and column 1.
    @inbounds for j in 1:self.buffer_pointer-1
        if (self.rowbuffer[j] > self.ndofs) || (self.rowbuffer[j] <= 0)
            self.rowbuffer[j] = 1
            self.matbuffer[j] = 0.0
        end
        if (self.colbuffer[j] > self.ndofs) || (self.colbuffer[j] <= 0)
            self.colbuffer[j] = 1
            self.matbuffer[j] = 0.0
        end
    end
    # We are making sure it is truly symmetric
    # @info "New"
    # r = view(self.rowbuffer, 1:self.buffer_pointer-1)
    # c = view(self.colbuffer, 1:self.buffer_pointer-1)
    # v = view(self.matbuffer, 1:self.buffer_pointer-1)
    # S = sparse(vcat(r, c), vcat(c, r), vcat(v, v), self.ndofs, self.ndofs); 
    # @info "Old"
    S = sparse(
        self.rowbuffer[1:self.buffer_pointer-1],
        self.colbuffer[1:self.buffer_pointer-1],
        self.matbuffer[1:self.buffer_pointer-1],
        self.ndofs,
        self.ndofs,
    )
    #  Now we need to construct the other triangle of the matrix. The diagonal
    #  will be duplicated.
    S = S + transpose(S)
    @inbounds for j in axes(S, 1)
        S[j, j] *= 0.5      # the diagonal is there twice; fix it;
    end
    self = SysmatAssemblerSparseSymm(zero(eltype(self.matbuffer))) # get rid of the buffers
    return S
end

"""
    SysmatAssemblerSparseDiag{T<:Number} <: AbstractSysmatAssembler

Assembler for a **symmetric square diagonal** matrix  assembled from symmetric
square diagonal matrices.

Warning: off-diagonal elements of the elementwise matrices will be ignored
during assembly!
"""
mutable struct SysmatAssemblerSparseDiag{IT, MBT, IBT} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length::IT
    matbuffer::MBT
    rowbuffer::IBT
    colbuffer::IBT
    buffer_pointer::IT
    ndofs::IT
end

"""
    SysmatAssemblerSparseDiag(zero::T=0.0) where {T<:Number}

Construct blank system matrix assembler for square diagonal matrices. The matrix
entries are of type `T`.

"""
function SysmatAssemblerSparseDiag(z= zero(FFlt))
    return SysmatAssemblerSparseDiag(0, FFlt[z], FInt[0], FInt[0], 0, 0)
end


"""
    startassembly!(self::SysmatAssemblerSparseDiag{T},
      elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt,
      ndofs::FInt, ignore2::FInt) where {T<:Number}

Start the assembly of a symmetric square global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.
- `elem_mat_nrows`= number of rows in typical element matrix,
- `elem_mat_ncols`= number of columns in a typical element matrix,
- `elem_mat_nmatrices`= number of element matrices,
- `ndofs_row`= Total number of equations in the row direction,
- `ndofs_col`= Total number of equations in the column direction.

The values stored in the buffers are initially undefined!

# Returns
- `self`: the modified assembler.
"""
function startassembly!(
    self::SysmatAssemblerSparseDiag,
    elem_mat_dim,
    ignore1,
    elem_mat_nmatrices,
    ndofs,
    ignore2,
)
    # Only reset the buffers if the pointer is 0.10
    if self.buffer_pointer < 1
        self.buffer_length = elem_mat_nmatrices * elem_mat_dim + 1
        resize!(self.rowbuffer, self.buffer_length)
        resize!(self.colbuffer, self.buffer_length)
        resize!(self.matbuffer, self.buffer_length)
        self.rowbuffer .= 1
        self.colbuffer .= 1
        self.matbuffer .= 0
        self.buffer_pointer = 1
        self.ndofs = ndofs
    end
    return self
end

"""
    assemble!(
        self::SysmatAssemblerSparseDiag,
        mat::MT,
        dofnums::IT,
        ignore
    ) where {MT, IT}

Assemble a square symmetric diagonal matrix.

- `dofnums` = the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are assumed to be the same).
- `mat` = **diagonal** square matrix
"""
function assemble!(
    self::SysmatAssemblerSparseDiag,
    mat::MT,
    dofnums::IT,
    ignore
) where {MT, IT}
    # Assembly of a square symmetric matrix.
    # The method assembles the lower triangle of the square symmetric matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums)
    ncolumns = nrows
    p = self.buffer_pointer
    @assert p + ncolumns <= self.buffer_length + 1
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns
        self.matbuffer[p] = mat[j, j] # serialized matrix
        self.rowbuffer[p] = dofnums[j]
        self.colbuffer[p] = dofnums[j]
        p = p + 1
    end
    self.buffer_pointer = p
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparseDiag)

Make a sparse symmetric square diagonal matrix.
"""
function makematrix!(self::SysmatAssemblerSparseDiag)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    @assert length(self.rowbuffer) >= self.buffer_pointer - 1
    @assert length(self.colbuffer) >= self.buffer_pointer - 1
    # Here we will go through the rows and columns, and whenever the row or
    # the column refer to indexes outside of the limits of the matrix, the
    # corresponding value will be set to 0 and assembled to row and column 1.
    @inbounds for j in 1:self.buffer_pointer-1
        if (self.rowbuffer[j] > self.ndofs) || (self.rowbuffer[j] <= 0)
            self.rowbuffer[j] = 1
            self.matbuffer[j] = 0.0
        end
        if (self.colbuffer[j] > self.ndofs) || (self.colbuffer[j] <= 0)
            self.colbuffer[j] = 1
            self.matbuffer[j] = 0.0
        end
    end
    S = sparse(
        self.rowbuffer[1:self.buffer_pointer-1],
        self.colbuffer[1:self.buffer_pointer-1],
        self.matbuffer[1:self.buffer_pointer-1],
        self.ndofs,
        self.ndofs,
    )
    self = SysmatAssemblerSparseDiag(zero(eltype(self.matbuffer))) # get rid of the buffers
    return S
end


# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================

"""
    AbstractSysvecAssembler

Abstract type of system vector assembler.
"""
abstract type AbstractSysvecAssembler end

"""
    startassembly!(self::SysvecAssembler{T}, ndofs_row::FInt) where {T<:Number}

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.

- `elem_mat_nmatrices` = number of element matrices expected to be processed
  during the assembly.
- `ndofs_row`= Total number of degrees of freedom.

# Returns
- `self`: the modified assembler.
"""
function startassembly!(self::SV, ndofs_row) where {SV<:AbstractSysvecAssembler} end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(
    self::SV,
    vec::MV,
    dofnums::IV,
) where {SV<:AbstractSysvecAssembler,MV,IV} end

"""
    makevector!(self::SysvecAssembler)

Make the global vector.
"""
function makevector!(self::SV) where {SV<:AbstractSysvecAssembler} end

"""
    SysvecAssembler

Assembler for the system vector.
"""
mutable struct SysvecAssembler{VBT, IT} <: AbstractSysvecAssembler
    F_buffer::VBT
    ndofs::IT
end

"""
    SysvecAssembler(z= zero(FFlt))

Construct blank system vector assembler. The vector entries are of type `T`
determined by the zero value.
"""
function SysvecAssembler(z= zero(FFlt))
    return SysvecAssembler([z], 1)
end

"""
    startassembly!(self::SysvecAssembler, ndofs_row::FInt) where {T<:Number}

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.

`ndofs_row`= Total number of degrees of freedom.

# Returns
- `self`: the modified assembler.
"""
function startassembly!(self::SysvecAssembler, ndofs_row)
    self.ndofs = ndofs_row
    resize!(self.F_buffer, self.ndofs)
    self.F_buffer .= 0
    return self
end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(
    self::SysvecAssembler,
    vec::MV,
    dofnums::IV,
) where {MV, IV}
    for i in eachindex(dofnums)
        gi = dofnums[i]
        if (0 < gi <= self.ndofs)
            self.F_buffer[gi] += vec[i]
        end
    end
end

"""
    makevector!(self::SysvecAssembler)

Make the global vector.
"""
function makevector!(self::SysvecAssembler)
    return deepcopy(self.F_buffer)
end


"""
    SysmatAssemblerSparseHRZLumpingSymm{IT, MBT, IBT} <: AbstractSysmatAssembler

Assembler for a **symmetric lumped square** matrix  assembled from  **symmetric square**
matrices. 

Reference: A note on mass lumping and related processes in the finite element method,
E. Hinton, T. Rock, O. C. Zienkiewicz, Earthquake Engineering & Structural Dynamics,
volume 4, number 3, 245--249, 1976.


!!! note 
    
    This assembler can compute and assemble diagonalized mass matrices.
    However, if the meaning of the entries of the mass matrix  differs
    (translation versus rotation), the mass matrices will not be computed
    correctly. Put bluntly: it can only be used for homogeneous mass matrices,
    all translation degrees of freedom, for instance. 
"""
mutable struct SysmatAssemblerSparseHRZLumpingSymm{IT, MBT, IBT} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length::IT
    matbuffer::MBT
    rowbuffer::IBT
    colbuffer::IBT
    buffer_pointer::IT
    ndofs::IT
end

"""
    SysmatAssemblerSparseHRZLumpingSymm(zer::T=0.0) where {T<:Number}

Construct blank system matrix assembler. The matrix entries are of type `T`.

"""
function SysmatAssemblerSparseHRZLumpingSymm(z= zero(FFlt))
    return SysmatAssemblerSparseHRZLumpingSymm(0, [z], Int[0], Int[0], 0, 0)
end

"""
    startassembly!(self::SysmatAssemblerSparseHRZLumpingSymm{T},
      elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt,
      ndofs::FInt, ignore2::FInt) where {T<:Number}

Start the assembly of a symmetric lumped square global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.
- `elem_mat_nrows`= number of rows in typical element matrix,
- `ignore1`= number of columns in a typical element matrix: equal to `elem_mat_nrows`,
- `elem_mat_nmatrices`= number of element matrices,
- `ndofs_row`= Total number of equations in the row direction,
- `ignore2`= Total number of equations in the column direction: equal to `ndofs_row`.

The values stored in the buffers are initially undefined!

# Returns
- `self`: the modified assembler.
"""
function startassembly!(
    self::SysmatAssemblerSparseHRZLumpingSymm,
    elem_mat_dim,
    ignore1,
    elem_mat_nmatrices,
    ndofs,
    ignore2,
)
    # Only reset the buffers if the pointer is 0.10
    if self.buffer_pointer < 1
        self.buffer_length = elem_mat_nmatrices * elem_mat_dim^2
        resize!(self.rowbuffer, self.buffer_length)
        resize!(self.colbuffer, self.buffer_length)
        resize!(self.matbuffer, self.buffer_length)
        self.rowbuffer .= 1
        self.colbuffer .= 1
        self.matbuffer .= 0
        self.buffer_pointer = 1
        self.ndofs = ndofs
    end
    return self
end

"""
    assemble!(
        self::SysmatAssemblerSparseHRZLumpingSymm,
        mat::MT,
        dofnums::IV,
        ignore::IV,
    ) where {MT, IV}

Assemble a HRZ-lumped square symmetric matrix.

Assembly of a HRZ-lumped square symmetric matrix. The method assembles the
scaled diagonal of the square symmetric matrix using the two vectors of
equation numbers for the rows and columns.
"""
function assemble!(
    self::SysmatAssemblerSparseHRZLumpingSymm,
    mat::MT,
    dofnums::IV,
    ignore::IV,
) where {MT, IV}
    nrows = length(dofnums)
    ncolumns = nrows
    p = self.buffer_pointer
    @assert p + ncolumns * nrows <= self.buffer_length + 1
    @assert size(mat) == (nrows, ncolumns)
    @assert nrows == ncolumns # only square matrices allowed
    # Now comes the lumping procedure
    em2 = sum(sum(mat, dims = 1)) # total mass times the number of space dimensions
    dem2 = zero(eltype(mat)) # total mass on the diagonal times the number of space dimensions
    for i in 1:nrows
        dem2 += mat[i, i]
    end
    ffactor = em2 / dem2 # total-element-mass compensation factor
    @inbounds for j in 1:ncolumns,  i in j:nrows
        if i == j # only the diagonal elements are assembled
            self.matbuffer[p] = mat[i, j] * ffactor # serialized matrix
            self.rowbuffer[p] = dofnums[i]
            self.colbuffer[p] = dofnums[j]
            p = p + 1
        end
    end
    self.buffer_pointer = p
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseHRZLumpingSymm{T}, mat::FMat{T},
        dofnums::FIntMat, ignore::FIntMat) where {T<:Number}

Assemble an HRZ-lumped square symmetric matrix.
"""
function assemble!(
    self::SysmatAssemblerSparseHRZLumpingSymm,
    mat::FMat{T},
    dofnums::FIntMat,
    ignore::FIntMat,
) where {T<:Number}
    return assemble!(self, mat, vec(dofnums), vec(ignore))
end

"""
    makematrix!(self::SysmatAssemblerSparseHRZLumpingSymm)

Make a sparse HRZ-lumped **symmetric square**  matrix.
"""
function makematrix!(self::SysmatAssemblerSparseHRZLumpingSymm)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    @assert length(self.rowbuffer) >= self.buffer_pointer - 1
    @assert length(self.colbuffer) >= self.buffer_pointer - 1
    # Here we will go through the rows and columns, and whenever the row or
    # the column refer to indexes outside of the limits of the matrix, the
    # corresponding value will be set to 0 and assembled to row and column 1.
    @inbounds for j in 1:self.buffer_pointer-1
        if (self.rowbuffer[j] > self.ndofs) || (self.rowbuffer[j] <= 0)
            self.rowbuffer[j] = 1
            self.matbuffer[j] = 0.0
        end
        if (self.colbuffer[j] > self.ndofs) || (self.colbuffer[j] <= 0)
            self.colbuffer[j] = 1
            self.matbuffer[j] = 0.0
        end
    end
    S = sparse(
        self.rowbuffer[1:self.buffer_pointer-1],
        self.colbuffer[1:self.buffer_pointer-1],
        self.matbuffer[1:self.buffer_pointer-1],
        self.ndofs,
        self.ndofs,
    )
    # Construct the other half of the matrix. (Even though this one really should be diagonal!)
    S = S + transpose(S)    # construct the other triangle
    @inbounds for j in axes(S, 1)
        S[j, j] /= 2.0      # the diagonal is there twice; fix it;
    end
    self = SysmatAssemblerSparse(0.0 * self.matbuffer[1])# get rid of the buffers
    return S
end


"""
    SysmatAssemblerReduced{T<:Number} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices.

!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix`.
"""
mutable struct SysmatAssemblerReduced{MT, TMT, IT} <: AbstractSysmatAssembler
    m::MT # reduced system matrix
    ndofs_row::IT
    ndofs_col::IT
    red_ndofs_row::IT
    red_ndofs_col::IT
    t::Matrix{TMT} # transformation matrix
end

function SysmatAssemblerReduced(t::TMT, z = zero(FFlt)) where {TMT}
    ndofs_row = ndofs_col = size(t, 1)
    red_ndofs_row = red_ndofs_col = size(t, 2)
    m = fill(z, red_ndofs_row, red_ndofs_col)
    return SysmatAssemblerReduced(
        m,
        ndofs_row,
        ndofs_col,
        red_ndofs_row,
        red_ndofs_col,
        t,
    )
end


function startassembly!(
    self::SysmatAssemblerReduced,
    elem_mat_nrows,
    elem_mat_ncols,
    elem_mat_nmatrices,
    ndofs_row,
    ndofs_col,
)
    @assert self.ndofs_row == ndofs_row
    @assert self.ndofs_col == ndofs_col
    self.m .= zero(eltype(self.m))
    return self
end

function assemble!(
    self::SysmatAssemblerReduced{T},
    mat::FMat{T},
    dofnums_row::Union{FIntVec,FIntMat},
    dofnums_col::Union{FIntVec,FIntMat},
) where {T<:Number}
    R, C = size(mat)
    @assert R == C "The matrix must be square"
    @assert dofnums_row == dofnums_col "The degree of freedom numbers must be the same for rows and columns"
    for i in 1:R
        gi = dofnums_row[i]
        if gi < 1
            mat[i, :] .= zero(T)
            dofnums_row[i] = 1
        end
    end
    for i in 1:C
        gi = dofnums_col[i]
        if gi < 1
            mat[:, i] .= zero(T)
            dofnums_col[i] = 1
        end
    end
    lt = self.t[dofnums_row, :]
    self.m .+= lt' * mat * lt
    return self
end

"""
    makematrix!(self::SysmatAssemblerReduced)

Make a sparse matrix.
"""
function makematrix!(self::SysmatAssemblerReduced)
    return self.m
end


end # end of module
