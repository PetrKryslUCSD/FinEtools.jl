"""
    AssemblyModule

Module for assemblers  of system matrices and vectors.
"""
module AssemblyModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import SparseArrays: sparse
import LinearAlgebra: diag

"""
    AbstractSysmatAssembler

Abstract type of system-matrix assembler.
"""
abstract type AbstractSysmatAssembler end;

"""
    SysmatAssemblerSparse{T<:Number} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices.
"""
mutable struct SysmatAssemblerSparse{T<:Number} <: AbstractSysmatAssembler
    buffer_length::FInt;
    matbuffer::Vector{T};
    rowbuffer::Vector{FInt};
    colbuffer::Vector{FInt};
    buffer_pointer::FInt;
    ndofs_row::FInt; ndofs_col::FInt;
end

"""
    SysmatAssemblerSparse(zero::T=0.0) where {T<:Number}

Construct blank system matrix assembler. The matrix entries are of type `T`.

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
"""
function SysmatAssemblerSparse(zero::T=0.0) where {T<:Number}
    return SysmatAssemblerSparse{T}(0,[zero],[0],[0],0,0,0)
end

"""
    startassembly!(self::SysmatAssemblerSparse{T},
      elem_mat_nrows::FInt, elem_mat_ncols::FInt, elem_mat_nmatrices::FInt,
      ndofs_row::FInt, ndofs_col::FInt) where {T<:Number}

Start the assembly of a global matrix.
"""
function startassembly!(self::SysmatAssemblerSparse{T},
  elem_mat_nrows::FInt, elem_mat_ncols::FInt, elem_mat_nmatrices::FInt,
  ndofs_row::FInt, ndofs_col::FInt) where {T<:Number}
    # Start the assembly of a global matrix.
    # The method makes buffers for matrix assembly. It must be called before
    # the first call to the method assemble.
    # elem_mat_nrows= number of rows in typical element matrix,
    # elem_mat_ncols= number of columns in a typical element matrix,
    # elem_mat_nmatrices= number of element matrices,
    # ndofs_row= Total number of equations in the row direction,
    # ndofs_col= Total number of equations in the column direction.
    self.buffer_length = elem_mat_nmatrices*elem_mat_nrows*elem_mat_ncols;
    self.rowbuffer = zeros(FInt,self.buffer_length);
    self.colbuffer = zeros(FInt,self.buffer_length);
    self.matbuffer = zeros(T,self.buffer_length);
    self.buffer_pointer = 1;
    self.ndofs_row = ndofs_row;
    self.ndofs_col = ndofs_col;
    return self
end

"""
    assemble!(self::SysmatAssemblerSparse{T}, mat::FMat{T},
      dofnums_row::FIntMat, dofnums_col::FIntMat) where {T<:Number}

Assemble a rectangular matrix.
"""
function assemble!(self::SysmatAssemblerSparse{T}, mat::FMat{T}, dofnums_row::FIntVec, dofnums_col::FIntVec) where {T<:Number}
    # Assembly of a rectangular matrix.
    # The method assembles a rectangular matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows=length(dofnums_row); ncolumns=length(dofnums_col);
    p = self.buffer_pointer
    @assert p+ncolumns*nrows <= self.buffer_length+1
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j=1:ncolumns
        @inbounds for i=1:nrows
            self.matbuffer[p] = mat[i,j] # serialized matrix
            self.rowbuffer[p] = dofnums_row[i];
            self.colbuffer[p] = dofnums_col[j];
            p=p+1
        end
    end
    self.buffer_pointer=p;
    return self
end

function assemble!(self::SysmatAssemblerSparse{T}, mat::FMat{T}, dofnums_row::FIntMat, dofnums_col::FIntMat) where {T<:Number}
    return assemble!(self, mat, vec(dofnums_row), vec(dofnums_col))
end

"""
    makematrix!(self::SysmatAssemblerSparse)

Make a sparse matrix.
"""
function makematrix!(self::SysmatAssemblerSparse)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    @assert length(self.rowbuffer) >= self.buffer_pointer-1
    @assert length(self.colbuffer) >= self.buffer_pointer-1
    @inbounds for j=1:self.buffer_pointer-1
        if (self.rowbuffer[j] > self.ndofs_row) || (self.rowbuffer[j] <= 0)
            self.rowbuffer[j]=self.ndofs_row+1;
        end
        if (self.colbuffer[j] > self.ndofs_col) || (self.colbuffer[j] <= 0)
            self.colbuffer[j]=self.ndofs_col+1;
        end
    end
    S = sparse(self.rowbuffer[1:self.buffer_pointer-1],
               self.colbuffer[1:self.buffer_pointer-1],
               self.matbuffer[1:self.buffer_pointer-1],
               self.ndofs_row+1, self.ndofs_col+1);
    self=SysmatAssemblerSparse(0.0*self.matbuffer[1])# get rid of the buffers
    return S[1:end-1,1:end-1]
end

"""
    SysmatAssemblerSparseSymm{T<:Number} <: AbstractSysmatAssembler

Assembler for a **symmetric square** matrix  assembled from symmetric square
matrices.
"""
mutable struct SysmatAssemblerSparseSymm{T<:Number} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length:: FInt;
    matbuffer::Vector{T};
    rowbuffer::Vector{FInt};
    colbuffer::Vector{FInt};
    buffer_pointer:: FInt;
    ndofs:: FInt;
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
function SysmatAssemblerSparseSymm(zero::T=0.0) where {T<:Number}
    return SysmatAssemblerSparseSymm{T}(0,[zero],[0],[0],0,0)
end

"""
    startassembly!(self::SysmatAssemblerSparseSymm{T},
      elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt,
      ndofs::FInt, ignore2::FInt) where {T<:Number}

Start the assembly of a symmetric square global matrix.
"""
function startassembly!(self::SysmatAssemblerSparseSymm{T},
  elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt,
  ndofs::FInt, ignore2::FInt) where {T<:Number}
    # Start the assembly of a global matrix.
    # The method makes buffers for matrix assembly. It must be called before
    # the first call to the method assemble.
    # elem_mat_nrows= number of rows in typical element matrix,
    # elem_mat_ncols= number of columns in a typical element matrix,
    # elem_mat_nmatrices= number of element matrices,
    # ndofs_row= Total number of equations in the row direction,
    # ndofs_col= Total number of equations in the column direction.
    self.buffer_length = elem_mat_nmatrices*elem_mat_dim^2;
    self.rowbuffer = zeros(FInt,self.buffer_length);
    self.colbuffer = zeros(FInt,self.buffer_length);
    self.matbuffer = zeros(T,self.buffer_length);
    self.buffer_pointer = 1;
    self.ndofs = ndofs;
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseSymm{T}, mat::FMat{T},  dofnums::FIntVec, ignore::FIntVec) where {T<:Number}

Assemble a square symmetric matrix.

`dofnums` are the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are the same).
"""
function assemble!(self::SysmatAssemblerSparseSymm{T}, mat::FMat{T},  dofnums::FIntVec, ignore::FIntVec) where {T<:Number}
    # Assembly of a square symmetric matrix.
    # The method assembles the lower triangle of the square symmetric matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows=length(dofnums); ncolumns=nrows;
    p = self.buffer_pointer
    @assert p+ncolumns*nrows <= self.buffer_length+1
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j=1:ncolumns
        @inbounds for i=j:nrows
            self.matbuffer[p] = mat[i,j] # serialized matrix
            self.rowbuffer[p] = dofnums[i];
            self.colbuffer[p] = dofnums[j];
            p=p+1
        end
    end
    self.buffer_pointer=p;
    return self
end

"""
	assemble!(self::SysmatAssemblerSparseSymm{T}, mat::FMat{T},
	  dofnums::FIntMat, ignore::FIntMat) where {T<:Number}

Assemble a square symmetric matrix.

`dofnums` are the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are the same).
"""
function assemble!(self::SysmatAssemblerSparseSymm{T}, mat::FMat{T}, dofnums::FIntMat, ignore::FIntMat) where {T<:Number}
    return assemble!(self, mat, vec(dofnums), vec(ignore))
end

"""
    makematrix!(self::SysmatAssemblerSparseSymm)

Make a sparse symmetric square matrix.
"""
function makematrix!(self::SysmatAssemblerSparseSymm)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    @assert length(self.rowbuffer) >= self.buffer_pointer-1
    @assert length(self.colbuffer) >= self.buffer_pointer-1
    @inbounds for j=1:self.buffer_pointer-1
        if (self.rowbuffer[j] > self.ndofs) || (self.rowbuffer[j] <= 0)
            self.rowbuffer[j]=self.ndofs+1;
        end
        if (self.colbuffer[j] > self.ndofs) || (self.colbuffer[j] <= 0)
            self.colbuffer[j]=self.ndofs+1;
        end
    end
    S = sparse(self.rowbuffer[1:self.buffer_pointer-1],
               self.colbuffer[1:self.buffer_pointer-1],
               self.matbuffer[1:self.buffer_pointer-1],
               self.ndofs+1, self.ndofs+1);   
    S = S+copy(transpose(S)); # This is to address the lack of a function for adding together a sparse matrix with a transpose of the sparse matrix; at the moment (January 2018), this defaults to the addition of two dense matrices and the result is DENSE. Hence a copy of the transpose needs to be made.
    #   S = S+transpose(S);    # construct the other triangle
    @inbounds for j=1:size(S,1)
        S[j,j]=S[j,j]/2.0;      # the diagonal is there twice; fix it;
    end
    self=SysmatAssemblerSparse(0.0*self.matbuffer[1])# get rid of the buffers
    S = S[1:end-1,1:end-1]
    return S
end


"""
    AbstractSysvecAssembler

Abstract type of system vector assembler.
"""
abstract type AbstractSysvecAssembler end;

"""
    startassembly!(self::SysvecAssembler{T},
      ndofs_row::FInt) where {T<:Number}

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.

`ndofs_row`= Total number of degrees of freedom.
"""
function startassembly!(self::SV,  ndofs_row::FInt) where {SV<:AbstractSysvecAssembler, T<:Number}
end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(self::SV, vec::MV, dofnums::D) where {SV<:AbstractSysvecAssembler, T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}
end

"""
    makevector!(self::SysvecAssembler)

Make the global vector.
"""
function makevector!(self::SV) where {SV<:AbstractSysvecAssembler}
end

"""
    SysvecAssembler

Assembler for the system vector.
"""
mutable struct SysvecAssembler{T<:Number} <: AbstractSysvecAssembler
    F_buffer::Vector{T};
    ndofs::FInt
end

"""
    SysvecAssembler(zero::T=0.0) where {T<:Number}

Construct blank system vector assembler. The vector entries are of type `T`.
"""
function SysvecAssembler(zero::T=0.0) where {T<:Number}
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
function startassembly!(self::SysvecAssembler{T},  ndofs_row::FInt) where {T<:Number}
    self.ndofs = ndofs_row
    self.F_buffer = zeros(T,self.ndofs);
end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(self::SysvecAssembler{T}, vec::MV,
  dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}
    for i = 1:length(dofnums)
        gi = dofnums[i];
        if (0 < gi <= self.ndofs)
            self.F_buffer[gi] = self.F_buffer[gi] + vec[i];
        end
    end
end

"""
    makevector!(self::SysvecAssembler)

Make the global vector.
"""
function makevector!(self::SysvecAssembler)
  return deepcopy(self.F_buffer);
end


"""
    SysmatAssemblerSparseHRZLumpingSymm{T<:Number} <: AbstractSysmatAssembler

Assembler for a **symmetric lumped square** matrix  assembled from  **symmetric square**
matrices. 

Reference: A note on mass lumping and related processes in the finite element method,
E. Hinton, T. Rock, O. C. Zienkiewicz, Earthquake Engineering & Structural Dynamics,
volume 4, number 3, 245--249, 1976.
}

"""
mutable struct SysmatAssemblerSparseHRZLumpingSymm{T<:Number} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length:: FInt;
    matbuffer::Vector{T};
    rowbuffer::Vector{FInt};
    colbuffer::Vector{FInt};
    buffer_pointer:: FInt;
    ndofs:: FInt;
end

"""
    SysmatAssemblerSparseHRZLumpingSymm(zer::T=0.0) where {T<:Number}

Construct blank system matrix assembler. The matrix entries are of type `T`.
"""
function SysmatAssemblerSparseHRZLumpingSymm(zer::T=0.0) where {T<:Number}
    return SysmatAssemblerSparseHRZLumpingSymm{T}(0,[zer],[0],[0],0,0)
end

"""
    startassembly!(self::SysmatAssemblerSparseHRZLumpingSymm{T},
      elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt,
      ndofs::FInt, ignore2::FInt) where {T<:Number}

Start the assembly of a symmetric lumped square global matrix.
"""
function startassembly!(self::SysmatAssemblerSparseHRZLumpingSymm{T}, elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt, ndofs::FInt,ignore2::FInt) where {T<:Number}
    # Start the assembly of a global matrix.
    # The method makes buffers for matrix assembly. It must be called before
    # the first call to the method assemble.
    # elem_mat_nrows= number of rows in typical element matrix,
    # elem_mat_ncols= number of columns in a typical element matrix,
    # elem_mat_nmatrices= number of element matrices,
    # ndofs_row= Total number of equations in the row direction,
    # ndofs_col= Total number of equations in the column direction.
    self.buffer_length = elem_mat_nmatrices*elem_mat_dim^2;
    self.rowbuffer = zeros(FInt,self.buffer_length);
    self.colbuffer = zeros(FInt,self.buffer_length);
    self.matbuffer = zeros(T,self.buffer_length);
    self.buffer_pointer = 1;
    self.ndofs = ndofs;
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseHRZLumpingSymm{T}, mat::FMat{T},
      dofnums::FIntMat, ignore::FIntMat) where {T<:Number}

Assemble a HRZ-lumped square symmetric matrix.
"""
function assemble!(self::SysmatAssemblerSparseHRZLumpingSymm{T}, mat::FMat{T},  dofnums::FIntVec, ignore::FIntVec) where {T<:Number}
    # Assembly of a HRZ-lumped square symmetric matrix.
    # The method assembles the scaled diagonal of the square symmetric matrix using the two vectors of equation numbers for the rows and columns.
    nrows=length(dofnums); ncolumns=nrows;
    p = self.buffer_pointer
    @assert p+ncolumns*nrows <= self.buffer_length+1
    @assert size(mat) == (nrows, ncolumns)
    # Now comes the lumping procedure
    em2 = sum(sum(mat, dims = 1));
    dem2 = sum(diag(mat));
    ffactor = em2/dem2 # total-element-mass compensation factor
    @inbounds for j=1:ncolumns
        @inbounds for i=j:nrows
            if i == j # only the diagonal elements are assembled
                self.matbuffer[p] = mat[i,j]*ffactor # serialized matrix
                self.rowbuffer[p] = dofnums[i];
                self.colbuffer[p] = dofnums[j];
                p=p+1
            end
        end
    end
    self.buffer_pointer=p;
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseHRZLumpingSymm{T}, mat::FMat{T},
        dofnums::FIntMat, ignore::FIntMat) where {T<:Number}

Assemble an HRZ-lumped square symmetric matrix.
"""
function assemble!(self::SysmatAssemblerSparseHRZLumpingSymm{T}, mat::FMat{T}, dofnums::FIntMat, ignore::FIntMat) where {T<:Number}
    return assemble!(self, mat, vec(dofnums), vec(ignore))
end

"""
    makematrix!(self::SysmatAssemblerSparseHRZLumpingSymm)

Make a sparse HRZ-lumped **symmetric square**  matrix.
"""
function makematrix!(self::SysmatAssemblerSparseHRZLumpingSymm)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    @assert length(self.rowbuffer) >= self.buffer_pointer-1
    @assert length(self.colbuffer) >= self.buffer_pointer-1
    @inbounds for j=1:self.buffer_pointer-1
        if (self.rowbuffer[j] > self.ndofs) || (self.rowbuffer[j] <= 0)
            self.rowbuffer[j]=self.ndofs+1;
        end
        if (self.colbuffer[j] > self.ndofs) || (self.colbuffer[j] <= 0)
            self.colbuffer[j]=self.ndofs+1;
        end
    end
    S = sparse(self.rowbuffer[1:self.buffer_pointer-1],
               self.colbuffer[1:self.buffer_pointer-1],
               self.matbuffer[1:self.buffer_pointer-1],
               self.ndofs+1, self.ndofs+1);   
    S = S+copy(transpose(S)); # This is to address the lack of a function for adding together a sparse matrix with a transpose of the sparse matrix; at the moment (January 2018), this defaults to the addition of two dense matrices and the result is DENSE. Hence a copy of the transpose needs to be made.
    #   S = S+transpose(S);    # construct the other triangle
    @inbounds for j=1:size(S,1)
        S[j,j]=S[j,j]/2.0;      # the diagonal is there twice; fix it;
    end
    self=SysmatAssemblerSparse(0.0*self.matbuffer[1])# get rid of the buffers
    S = S[1:end-1,1:end-1]
    return S
end


end
