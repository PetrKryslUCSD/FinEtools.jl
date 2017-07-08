"""
    AssemblyModule

Module for assemblers  of system matrices and vectors.  
"""
module AssemblyModule

export SysmatAssemblerBase, SysmatAssemblerSparse, SysmatAssemblerSparseSymm
export startassembly!, assemble!, makematrix!
export SysvecAssemblerBase, SysvecAssembler
export startassembly!, assemble!, makevector!

using FinEtools.FTypesModule

const  inv_dofnum=0;            # invalid degree of freedom number -- no equation

abstract type SysmatAssemblerBase end;


mutable struct SysmatAssemblerSparse{T<:Number} <: SysmatAssemblerBase
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length::FInt;
    matbuffer::FVec{T};
    rowbuffer::FIntVec;
    colbuffer::FIntVec;
    buffer_pointer::FInt;
    ndofs_row::FInt; ndofs_col::FInt;
end

# AssemblyModule.SysmatAssemblerSparse()
function SysmatAssemblerSparse{T<:Number}(zero::T=0.0)
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
    self.ndofs_row =ndofs_row;
    self.ndofs_col =ndofs_col;
    return self
end

"""
    assemble!(self::SysmatAssemblerSparse{T}, mat::FMat{T},
      dofnums_row::FIntMat, dofnums_col::FIntMat) where {T<:Number}

Assembly of a rectangular matrix.
"""
function assemble!(self::SysmatAssemblerSparse{T}, mat::FMat{T},
  dofnums_row::FIntMat, dofnums_col::FIntMat) where {T<:Number}
    # Assembly of a rectangular matrix.
    # The method assembles a rectangular matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows=length(dofnums_row); ncolumns=length(dofnums_col);
    p = self.buffer_pointer
    @assert p+ncolumns*nrows <= self.buffer_length+1
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

"""
    makematrix!(self::SysmatAssemblerSparse)

Make a sparse matrix.
"""
function makematrix!(self::SysmatAssemblerSparse)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    @inbounds for j=1:self.buffer_pointer-1
        if self.rowbuffer[j] == inv_dofnum
            self.rowbuffer[j]=self.ndofs_row+1;
        end
        if self.colbuffer[j] == inv_dofnum
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
    SysmatAssemblerSparseSymm{T<:Number} <: SysmatAssemblerBase

Assembler for a SYMMETRIC square matrix  assembled from  SYMMETRIC squared
matrices.
"""
mutable struct SysmatAssemblerSparseSymm{T<:Number} <: SysmatAssemblerBase
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length:: FInt;
    matbuffer::FVec{T};
    rowbuffer::FIntVec;
    colbuffer::FIntVec;
    buffer_pointer:: FInt;
    ndofs:: FInt;
end

# AssemblyModule.SysmatAssemblerSparse()
function SysmatAssemblerSparseSymm{T<:Number}(zero::T=0.0)
    return SysmatAssemblerSparseSymm{T}(0,[zero],[0],[0],0,0)
end

"""
    startassembly!(self::SysmatAssemblerSparseSymm{T},
      elem_mat_dim::FInt, ignore1::FInt, elem_mat_nmatrices::FInt,
      ndofs::FInt, ignore2::FInt) where {T<:Number}

Start the assembly of a SYMMETRIC square global matrix.
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
    self.ndofs =ndofs;
    return self
end

"""
    assemble!(self::SysmatAssemblerSparseSymm{T}, mat::FMat{T},
      dofnums::FIntMat, ignore::FIntMat) where {T<:Number}

Assembly of a square symmetric matrix.
"""
function assemble!(self::SysmatAssemblerSparseSymm{T}, mat::FMat{T},
  dofnums::FIntMat, ignore::FIntMat) where {T<:Number}
    # Assembly of a square symmetric matrix.
    # The method assembles the lower triangle of the square symmetric matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows=length(dofnums); ncolumns=nrows;
    p = self.buffer_pointer
    @assert p+ncolumns*nrows <= self.buffer_length+1
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
    makematrix!(self::SysmatAssemblerSparseSymm)

Make a sparse SYMMETRIC SQUARE matrix.
"""
function makematrix!(self::SysmatAssemblerSparseSymm)
    # Make a sparse matrix.
    # The method makes a sparse matrix from the assembly buffers.
    @inbounds for j=1:self.buffer_pointer-1
        if self.rowbuffer[j] == inv_dofnum
            self.rowbuffer[j]=self.ndofs+1;
        end
        if self.colbuffer[j] == inv_dofnum
            self.colbuffer[j]=self.ndofs+1;
        end
    end
    S = sparse(self.rowbuffer[1:self.buffer_pointer-1],
               self.colbuffer[1:self.buffer_pointer-1],
               self.matbuffer[1:self.buffer_pointer-1],
               self.ndofs+1, self.ndofs+1);
    S = S+S';    # construct the other triangle
    @inbounds for j=1:size(S,1)
        S[j,j]=S[j,j]/2.0;      # the diagonal is there twice; fix it;
    end
    self=SysmatAssemblerSparse(0.0*self.matbuffer[1])# get rid of the buffers
    return S[1:end-1,1:end-1]
end



abstract type SysvecAssemblerBase end;

mutable struct SysvecAssembler{T<:Number} <: SysvecAssemblerBase
    # % The class sysvec_assembler is for assembling of a system
    # % column vector from elementwise vectors.
    # 	% Fixed degrees of freedom numbers are given this value:
    #      % it indicates that this is not a valid free  degree of freedom number.
    F_buffer::FVec{T};
end

function SysvecAssembler{T<:Number}(zero::T=0.0)
    return SysvecAssembler([zero])
end

"""
    startassembly!(self::SysvecAssembler{T},
      ndofs_row::FInt) where {T<:Number}

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.
ndofs_row= Total number of degrees of freedom.
"""
function startassembly!(self::SysvecAssembler{T},
  ndofs_row::FInt) where {T<:Number}
    self.F_buffer= zeros(T,ndofs_row);
end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assembly of elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(self::SysvecAssembler{T}, vec::MV,
  dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}
    for i = 1:length(dofnums)
        gi = dofnums[i];
        if (gi != inv_dofnum)
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


end
