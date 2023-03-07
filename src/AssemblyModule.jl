"""
    AssemblyModule

Module for assemblers  of system matrices and vectors.
"""
module AssemblyModule

__precompile__(true)

using ..FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using SparseArrays: sparse, spzeros
using LinearAlgebra: diag

"""
    AbstractSysmatAssembler

Abstract type of system-matrix assembler.
"""
abstract type AbstractSysmatAssembler end

"""
    SysmatAssemblerSparse{T<:Number} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices.

!!! note

    All fields of the datatype are private. No need to access them directly.
"""
mutable struct SysmatAssemblerSparse{T<:Number} <: AbstractSysmatAssembler
    buffer_length::FInt
    matbuffer::Vector{T}
    rowbuffer::Vector{FInt}
    colbuffer::Vector{FInt}
    buffer_pointer::FInt
    ndofs_row::FInt
    ndofs_col::FInt
    nomatrixresult::Bool
end

"""
    SysmatAssemblerSparse(zero::T=0.0, nomatrixresult = false) where {T<:Number}

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
function SysmatAssemblerSparse(zero::T = 0.0, nomatrixresult = false) where {T<:Number}
    return SysmatAssemblerSparse{T}(0, [zero], [0], [0], 0, 0, 0, nomatrixresult)
end

"""
    startassembly!(self::SysmatAssemblerSparse{T},
      elem_mat_nrows::FInt, elem_mat_ncols::FInt, elem_mat_nmatrices::FInt,
      ndofs_row::FInt, ndofs_col::FInt) where {T<:Number}

Start the assembly of a global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.
- `elem_mat_nrows`= number of rows in typical element matrix,
- `elem_mat_ncols`= number of columns in a typical element matrix,
- `elem_mat_nmatrices`= number of element matrices,
- `ndofs_row`= Total number of equations in the row direction,
- `ndofs_col`= Total number of equations in the column direction.

The values stored in the buffers are initially undefined!
"""
function startassembly!(
    self::SysmatAssemblerSparse{T},
    elem_mat_nrows::FInt,
    elem_mat_ncols::FInt,
    elem_mat_nmatrices::FInt,
    ndofs_row::FInt,
    ndofs_col::FInt,
) where {T<:Number}
    self.buffer_length = elem_mat_nmatrices * elem_mat_nrows * elem_mat_ncols
    self.rowbuffer = Array{FInt,1}(undef, self.buffer_length)
    self.colbuffer = Array{FInt,1}(undef, self.buffer_length)
    self.matbuffer = Array{T,1}(undef, self.buffer_length)
    self.buffer_pointer = 1
    self.ndofs_row = ndofs_row
    self.ndofs_col = ndofs_col
    return self
end

"""
    assemble!(self::SysmatAssemblerSparse{T}, mat::FMat{T}, dofnums_row::Union{FIntVec, FIntMat}, dofnums_col::Union{FIntVec, FIntMat}) where {T<:Number}

Assemble a rectangular matrix.
"""
function assemble!(
    self::SysmatAssemblerSparse{T},
    mat::FMat{T},
    dofnums_row::Union{FIntVec,FIntMat},
    dofnums_col::Union{FIntVec,FIntMat},
) where {T<:Number}
    # Assembly of a rectangular matrix.
    # The method assembles a rectangular matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    p = self.buffer_pointer
    @assert p + ncolumns * nrows <= self.buffer_length + 1
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns
        @inbounds for i in 1:nrows
            self.matbuffer[p] = mat[i, j] # serialized matrix
            self.rowbuffer[p] = dofnums_row[i]
            self.colbuffer[p] = dofnums_col[j]
            p = p + 1
        end
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
    SysmatAssemblerSparseSymm{T<:Number} <: AbstractSysmatAssembler

Assembler for a **symmetric square** matrix  assembled from symmetric square
matrices.
"""
mutable struct SysmatAssemblerSparseSymm{T<:Number} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length::FInt
    matbuffer::Vector{T}
    rowbuffer::Vector{FInt}
    colbuffer::Vector{FInt}
    buffer_pointer::FInt
    ndofs::FInt
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
[]
"""
function SysmatAssemblerSparseSymm(zero::T = 0.0) where {T<:Number}
    return SysmatAssemblerSparseSymm{T}(0, [zero], [0], [0], 0, 0)
end

"""
    startassembly!(self::SysmatAssemblerSparseSymm{T},
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
"""
function startassembly!(
    self::SysmatAssemblerSparseSymm{T},
    elem_mat_dim::FInt,
    ignore1::FInt,
    elem_mat_nmatrices::FInt,
    ndofs::FInt,
    ignore2::FInt,
) where {T<:Number}
    self.buffer_length = elem_mat_nmatrices * elem_mat_dim^2
    self.rowbuffer = Array{FInt,1}(undef, self.buffer_length)
    self.colbuffer = Array{FInt,1}(undef, self.buffer_length)
    self.matbuffer = Array{T,1}(undef, self.buffer_length)
    self.buffer_pointer = 1
    self.ndofs = ndofs
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseSymm{T}, mat::FMat{T},  dofnums::FIntVec, ignore::FIntVec) where {T<:Number}

Assemble a square symmetric matrix.

`dofnums` are the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are assumed to be the same).
"""
function assemble!(
    self::SysmatAssemblerSparseSymm{T},
    mat::FMat{T},
    dofnums::Union{FIntVec,FIntMat},
    ignore::Union{FIntVec,FIntMat},
) where {T<:Number}
    # Assembly of a square symmetric matrix.
    # The method assembles the lower triangle of the square symmetric matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums)
    ncolumns = nrows
    p = self.buffer_pointer
    @assert p + ncolumns * nrows <= self.buffer_length + 1
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns
        @inbounds for i in j:nrows
            self.matbuffer[p] = mat[i, j] # serialized matrix
            self.rowbuffer[p] = dofnums[i]
            self.colbuffer[p] = dofnums[j]
            p = p + 1
        end
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
mutable struct SysmatAssemblerSparseDiag{T<:Number} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length::FInt
    matbuffer::Vector{T}
    rowbuffer::Vector{FInt}
    colbuffer::Vector{FInt}
    buffer_pointer::FInt
    ndofs::FInt
end

"""
    SysmatAssemblerSparseDiag(zero::T=0.0) where {T<:Number}

Construct blank system matrix assembler for symmetric matrices. The matrix
entries are of type `T`.

# Example

This is how a symmetric sparse matrix is assembled from two square dense matrices.
```
    a = SysmatAssemblerSparseDiag(0.0)                                                        
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
[]
"""
function SysmatAssemblerSparseDiag(zero::T = 0.0) where {T<:Number}
    return SysmatAssemblerSparseDiag{T}(0, [zero], [0], [0], 0, 0)
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
"""
function startassembly!(
    self::SysmatAssemblerSparseDiag{T},
    elem_mat_dim::FInt,
    ignore1::FInt,
    elem_mat_nmatrices::FInt,
    ndofs::FInt,
    ignore2::FInt,
) where {T<:Number}
    self.buffer_length = elem_mat_nmatrices * elem_mat_dim + 1
    self.rowbuffer = Array{FInt,1}(undef, self.buffer_length)
    self.colbuffer = Array{FInt,1}(undef, self.buffer_length)
    self.matbuffer = Array{T,1}(undef, self.buffer_length)
    self.buffer_pointer = 1
    self.ndofs = ndofs
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseDiag{T}, mat::FMat{T},  dofnums::FIntVec, ignore::FIntVec) where {T<:Number}

Assemble a square symmetric diagonal matrix.

- `dofnums` = the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are assumed to be the same).
- `mat` = diagonal square matrix
"""
function assemble!(
    self::SysmatAssemblerSparseDiag{T},
    mat::FMat{T},
    dofnums::Union{FIntVec,FIntMat},
    ignore::Union{FIntVec,FIntMat},
) where {T<:Number}
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

"""
    SysmatAssemblerReduced{T<:Number} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices.

!!! note

    All fields of the datatype are private. No need to access them directly.
"""
mutable struct SysmatAssemblerReduced{T<:Number} <: AbstractSysmatAssembler
    m::Matrix{T} # reduced system matrix
    ndofs_row::FInt
    ndofs_col::FInt
    red_ndofs_row::FInt
    red_ndofs_col::FInt
    t::Matrix{T} # transformation matrix
end

function SysmatAssemblerReduced(t::Matrix{T}) where {T<:Number}
    ndofs_row = ndofs_col = size(t, 1)
    red_ndofs_row = red_ndofs_col = size(t, 2)
    m = fill(zero(T), red_ndofs_row, red_ndofs_col)
    return SysmatAssemblerReduced{T}(
        m,
        ndofs_row,
        ndofs_col,
        red_ndofs_row,
        red_ndofs_col,
        t,
    )
end


function startassembly!(
    self::SysmatAssemblerReduced{T},
    elem_mat_nrows::FInt,
    elem_mat_ncols::FInt,
    elem_mat_nmatrices::FInt,
    ndofs_row::FInt,
    ndofs_col::FInt,
) where {T<:Number}
    @assert self.ndofs_row == ndofs_row
    @assert self.ndofs_col == ndofs_col
    self.m[:] .= zero(T)
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
    startassembly!(self::SysvecAssembler{T},
      ndofs_row::FInt) where {T<:Number}

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.

`ndofs_row`= Total number of degrees of freedom.
"""
function startassembly!(self::SV, ndofs_row::FInt) where {SV<:AbstractSysvecAssembler} end

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
    dofnums::D,
) where {SV<:AbstractSysvecAssembler,T<:Number,MV<:AbstractArray{T},D<:AbstractArray{FInt}} end

"""
    makevector!(self::SysvecAssembler)

Make the global vector.
"""
function makevector!(self::SV) where {SV<:AbstractSysvecAssembler} end

"""
    SysvecAssembler

Assembler for the system vector.
"""
mutable struct SysvecAssembler{T<:Number} <: AbstractSysvecAssembler
    F_buffer::Vector{T}
    ndofs::FInt
end

"""
    SysvecAssembler(zero::T=0.0) where {T<:Number}

Construct blank system vector assembler. The vector entries are of type `T`.
"""
function SysvecAssembler(zero::T = 0.0) where {T<:Number}
    return SysvecAssembler([zero], 1)
end

"""
    startassembly!(self::SysvecAssembler{T},
      ndofs_row::FInt) where {T<:Number}

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.

`ndofs_row`= Total number of degrees of freedom.
"""
function startassembly!(self::SysvecAssembler{T}, ndofs_row::FInt) where {T<:Number}
    self.ndofs = ndofs_row
    self.F_buffer = zeros(T, self.ndofs)
end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(
    self::SysvecAssembler{T},
    vec::MV,
    dofnums::D,
) where {T<:Number,MV<:AbstractArray{T},D<:AbstractArray{FInt}}
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
    SysmatAssemblerSparseHRZLumpingSymm{T<:Number} <: AbstractSysmatAssembler

Assembler for a **symmetric lumped square** matrix  assembled from  **symmetric square**
matrices. 

Reference: A note on mass lumping and related processes in the finite element method,
E. Hinton, T. Rock, O. C. Zienkiewicz, Earthquake Engineering & Structural Dynamics,
volume 4, number 3, 245--249, 1976.
}


!!! note 
    
    This assembler can compute and assemble diagonalized mass matrices.
    However, if the meaning of the entries of the mass matrix  differs
    (translation versus rotation), the mass matrices will not be computed
    correctly. Put bluntly: it can only be used for homogeneous mass matrices,
    all translation degrees of freedom, for instance. 
"""
mutable struct SysmatAssemblerSparseHRZLumpingSymm{T<:Number} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length::FInt
    matbuffer::Vector{T}
    rowbuffer::Vector{FInt}
    colbuffer::Vector{FInt}
    buffer_pointer::FInt
    ndofs::FInt
end

"""
    SysmatAssemblerSparseHRZLumpingSymm(zer::T=0.0) where {T<:Number}

Construct blank system matrix assembler. The matrix entries are of type `T`.

"""
function SysmatAssemblerSparseHRZLumpingSymm(zer::T = 0.0) where {T<:Number}
    return SysmatAssemblerSparseHRZLumpingSymm{T}(0, [zer], [0], [0], 0, 0)
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
"""
function startassembly!(
    self::SysmatAssemblerSparseHRZLumpingSymm{T},
    elem_mat_dim::FInt,
    ignore1::FInt,
    elem_mat_nmatrices::FInt,
    ndofs::FInt,
    ignore2::FInt,
) where {T<:Number}
    self.buffer_length = elem_mat_nmatrices * elem_mat_dim^2
    self.rowbuffer = Array{FInt,1}(undef, self.buffer_length)
    self.colbuffer = Array{FInt,1}(undef, self.buffer_length)
    self.matbuffer = Array{T,1}(undef, self.buffer_length)
    self.buffer_pointer = 1
    self.ndofs = ndofs
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseHRZLumpingSymm{T}, mat::FMat{T},
      dofnums::FIntMat, ignore::FIntMat) where {T<:Number}

Assemble a HRZ-lumped square symmetric matrix.

Assembly of a HRZ-lumped square symmetric matrix. The method assembles the
scaled diagonal of the square symmetric matrix using the two vectors of
equation numbers for the rows and columns.
"""
function assemble!(
    self::SysmatAssemblerSparseHRZLumpingSymm{T},
    mat::FMat{T},
    dofnums::FIntVec,
    ignore::FIntVec,
) where {T<:Number}

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
    @inbounds for j in 1:ncolumns
        @inbounds for i in j:nrows
            if i == j # only the diagonal elements are assembled
                self.matbuffer[p] = mat[i, j] * ffactor # serialized matrix
                self.rowbuffer[p] = dofnums[i]
                self.colbuffer[p] = dofnums[j]
                p = p + 1
            end
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
    self::SysmatAssemblerSparseHRZLumpingSymm{T},
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
    SysmatAssemblerSparseThr{T<:Number} <: AbstractSysmatAssembler

Type for multithreaded assembling of a sparse global matrix from elementwise matrices.

!!! note

    All fields of the datatype are private. No need to access them directly.
"""
mutable struct SysmatAssemblerSparseThr{T<:Number} <: AbstractSysmatAssembler
    matbuffer::Vector{Vector{T}}
    rowbuffer::Vector{Vector{FInt}}
    colbuffer::Vector{Vector{FInt}}
    thread_done::Vector{Bool}
    ndofs_row::FInt
    ndofs_col::FInt
    nomatrixresult::Bool
end

"""
    SysmatAssemblerSparseThr(zero::T=0.0, nomatrixresult = false) where {T<:Number}

Construct blank system matrix assembler.

This assembler is like `SysmatAssemblerSparse`, except it can be used with
multiple threads.
"""
function SysmatAssemblerSparseThr(zero::T = 0.0, nomatrixresult = false) where {T<:Number}
    return SysmatAssemblerSparseThr{T}(Vector{Vector{T}}[], Vector{Vector{FInt}}[], Vector{Vector{FInt}}[], Vector{Bool}[], 0, 0, nomatrixresult)
end

"""
    startassembly!(self::SysmatAssemblerSparseThr{T},
      elem_mat_nrows::FInt, elem_mat_ncols::FInt, elem_mat_nmatrices::FInt,
      ndofs_row::FInt, ndofs_col::FInt) where {T<:Number}

Start the assembly of a global matrix.

This assembler is like `SysmatAssemblerSparse`, except it can be used with
multiple threads.
"""
function startassembly!(
    self::SysmatAssemblerSparseThr{T},
    elem_mat_nrows::FInt,
    elem_mat_ncols::FInt,
    elem_mat_nmatrices::FInt,
    ndofs_row::FInt,
    ndofs_col::FInt,
) where {T<:Number}
    nth = Base.Threads.nthreads()
    n = Int(ceil(round(elem_mat_nmatrices * elem_mat_nrows * elem_mat_ncols / nth)))
    self.rowbuffer = Vector{Vector{FInt}}[]
    self.colbuffer = Vector{Vector{FInt}}[]
    self.matbuffer = Vector{Vector{T}}[]
    self.thread_done = fill(false, nth)
    for th in 1:nth
        push!(self.rowbuffer, Vector{FInt}[])
        sizehint!(self.rowbuffer[th], n)
        push!(self.colbuffer, Vector{FInt}[])
        sizehint!(self.colbuffer[th], n)
        push!(self.matbuffer, Vector{T}[])
        sizehint!(self.matbuffer[th], n)
    end
    self.ndofs_row = ndofs_row
    self.ndofs_col = ndofs_col
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseThr{T}, mat::FMat{T}, dofnums_row::Union{FIntVec, FIntMat}, dofnums_col::Union{FIntVec, FIntMat}) where {T<:Number}

Assemble a rectangular matrix.

This assembler is like `SysmatAssemblerSparse`, except it can be used with
multiple threads.
"""
function assemble!(
    self::SysmatAssemblerSparseThr{T},
    mat::FMat{T},
    dofnums_row::Union{FIntVec,FIntMat},
    dofnums_col::Union{FIntVec,FIntMat},
) where {T<:Number}
    nth = Base.Threads.nthreads()
    th = Base.Threads.threadid()
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    @inbounds for j in 1:ncolumns, i in 1:nrows
        push!(self.matbuffer[th], mat[i, j]) # serialized matrix
        push!(self.rowbuffer[th], dofnums_row[i])
        push!(self.colbuffer[th], dofnums_col[j])
    end
    @assert length(self.matbuffer[th]) == length(self.rowbuffer[th])
    @assert length(self.matbuffer[th]) == length(self.colbuffer[th])
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparseThr)

Make a sparse matrix.

This assembler is like `SysmatAssemblerSparse`, except it can be used with
multiple threads.
"""
function makematrix!(self::SysmatAssemblerSparseThr{T}) where {T<:Number}
    # Here we will go through the rows and columns, and whenever the row or
    # the column refer to indexes outside of the limits of the matrix, the
    # corresponding value will be set to 0 and assembled to row and column 1.
    nth = Base.Threads.nthreads()
    th = Base.Threads.threadid()
    @inbounds for j in 1:length(self.rowbuffer[th])
        if (self.rowbuffer[th][j] > self.ndofs_row) || (self.rowbuffer[th][j] <= 0)
            self.rowbuffer[th][j] = 1
            self.matbuffer[th][j] = 0.0
        end
        if (self.colbuffer[th][j] > self.ndofs_col) || (self.colbuffer[th][j] <= 0)
            self.colbuffer[th][j] = 1
            self.matbuffer[th][j] = 0.0
        end
    end
    self.thread_done[th] = true
    all_done = true
    for k in eachindex(self.thread_done)
        @assert length(self.matbuffer[k]) == length(self.rowbuffer[k])
        @assert length(self.matbuffer[k]) == length(self.colbuffer[k])
        all_done = all_done && (
            self.thread_done[k] || (
        isempty(self.rowbuffer[k]) && isempty(self.colbuffer[k]) && isempty(self.matbuffer[k])))
    end
    if self.nomatrixresult || (!all_done)
        # No actual sparse matrix is returned. The entire result of the assembly
        # is preserved in the assembler buffers.
        return spzeros(self.ndofs_row, self.ndofs_col)
    else
        # The sparse matrix is constructed and returned. The  buffers used for
        # the assembly are cleared.
        S = sparse(
            cat(self.rowbuffer..., dims=1),
            cat(self.colbuffer..., dims=1),
            cat(self.matbuffer..., dims=1),
            self.ndofs_row,
            self.ndofs_col,
        )
        self = SysmatAssemblerSparseThr(zero(T)) # get rid of the buffers
        return S
    end
end

end
