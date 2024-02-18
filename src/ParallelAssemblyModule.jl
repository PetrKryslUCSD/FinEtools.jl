module ParallelAssemblyModule

using ..FEMMBaseModule: AbstractFEMM, finite_elements, iselementbased
using ..FieldModule: ndofs, nalldofs
using ..NodalFieldModule: NodalField
using ..AssemblyModule:
    AbstractSysmatAssembler, startassembly!, makematrix!, setnomatrixresult, expectedntriples
using ..FESetModule: nodesperelem

function _update_buffer_range(buffer_length, iend)
    istart = iend + 1
    iend = iend + buffer_length
    buffer_range = istart:iend
    return buffer_range, iend
end

function _task_local_assembler(
    AssemblerType::Type{AT},
    buffer_range,
    matbuffer,
    rowbuffer,
    colbuffer,
    row_nalldofs,
    col_nalldofs,
) where {AT}
    buffer_length = maximum(buffer_range) - minimum(buffer_range) + 1
    buffer_pointer = 1
    nomatrixresult = true
    force_init = false
    return AssemblerType(
        buffer_length,
        view(matbuffer, buffer_range),
        view(rowbuffer, buffer_range),
        view(colbuffer, buffer_range),
        buffer_pointer,
        row_nalldofs,
        col_nalldofs,
        nomatrixresult,
        force_init,
    )
end

function _check_femm_compatibility(femms::AbstractArray{FEMM,1}) where {FEMM<:AbstractFEMM}
    for j in eachindex(femms)
        iselementbased(femms[j]) || error("FEMM is not element-based")
        (nameof(typeof(femms[j])) === nameof(typeof(femms[1]))) ||
            error("All FEMMs must be of the same type")
        (
            nameof(typeof(finite_elements(femms[j]))) ===
            nameof(typeof(finite_elements(femms[1])))
        ) || error("All finite elements must be of the same type")
    end
    return true
end

"""
    make_assembler(
        femms::AbstractArray{FEMM, 1},
        AssemblerType::Type{AT},
        uf::NodalField{UFT},
    ) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, UFT<:Number}
    
Make the global assembler.    

- `femms`: array of finite element machines, one for each subdomain
- `AssemblerType`: type of the global assembler, unparametrized
- `uf`: nodal field of the unknowns
"""
function make_assembler(
    femms::AbstractArray{FEMM,1},
    AssemblerType::Type{AT},
    uf::NodalField{UFT},
) where {FEMM<:AbstractFEMM,AT<:AbstractSysmatAssembler,UFT<:Number}
    _check_femm_compatibility(femms)
    elem_mat_nmatrices = count.(finite_elements.(femms))
    # Count the total number of COO triples
    nne = nodesperelem(finite_elements(femms[1]))
    ndn = ndofs(uf)
    elem_mat_nrows = nne * ndn
    elem_mat_ncols = nne * ndn
    ndofs_row = nalldofs(uf)
    ndofs_col = nalldofs(uf)
    # Start the assembly
    assembler = startassembly!(
        AssemblerType(zero(UFT)),
        elem_mat_nrows,
        elem_mat_ncols,
        sum(elem_mat_nmatrices),
        ndofs_row,
        ndofs_col,
    )
    setnomatrixresult(assembler, true)
    return assembler
end

"""
    start_assembler!(
        assembler::AT
    ) where {AT<:AbstractSysmatAssembler}

Start parallel assembly on the global assembler.

The assembler is informed that no matrix result should be computed. All the data
should be preserved in the assembly buffers.
"""
function start_assembler!(assembler::AT) where {AT<:AbstractSysmatAssembler}
    assembler = startassembly!(
        assembler,
        1,
        1,
        assembler._buffer_length,
        assembler._row_nalldofs,
        assembler._col_nalldofs,
    )
    setnomatrixresult(assembler, true)
    return assembler
end


"""
    make_matrix!(
        assembler::AT
    ) where {AT<:AbstractSysmatAssembler}

Make the global matrix, using the data collected in the assembler.
"""
function make_matrix!(assembler::AT) where {AT<:AbstractSysmatAssembler}
    # At this point all data is in the buffer
    assembler._buffer_pointer = assembler._buffer_length + 1
    setnomatrixresult(assembler, false)
    return makematrix!(assembler)
end

"""
    make_task_assemblers(
        femms::AbstractArray{FEMM, 1},
        assembler::AT,
        AssemblerType::Type{AT},
        uf::NodalField{UFT}
    ) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, UFT<:Number}

Make task assemblers for the subdomains.

- `femms`: array of finite element machines, one for each subdomain
- `assembler`: the global assembler 
- `AssemblerType`: type of the global assembler, unparametrized
- `uf`: nodal field of the unknowns

Return an array of the assemblers. Each of the assemblers only refers to a view
of the buffers of the global assembler, hence they are reasonably cheap to
construct.
"""
function make_task_assemblers(
    femms::AbstractArray{FEMM,1},
    assembler::AT,
    AssemblerType::Type{AT},
    uf::NodalField{UFT},
) where {FEMM<:AbstractFEMM,AT<:AbstractSysmatAssembler,UFT<:Number}
    _check_femm_compatibility(femms)
    # The number of tasks is equal to the number of the FEMMs
    ntasks = length(femms)
    ntasks >= 1 || error("Number of tasks must be >= 1")
    elem_mat_nmatrices = count.(finite_elements.(femms))
    # Check compatibility of the FEMMs: they need to be the same "type" -- all
    # conductivity, for instance, and all using the same type of finite element 
    elem_mat_nmatrices = count.(finite_elements.(femms))
    # Count the total number of COO triples
    nne = nodesperelem(finite_elements(femms[1]))
    ndn = ndofs(uf)
    elem_mat_nrows = nne * ndn
    elem_mat_ncols = nne * ndn
    triples = [expectedntriples(assembler, elem_mat_nrows, elem_mat_ncols, n) for n in elem_mat_nmatrices]
    # Make there task assemblers
    assemblers = AssemblerType[]
    iend = 0
    for j in eachindex(femms)
        buffer_range, iend = _update_buffer_range(triples[j], iend)
        push!(
            assemblers,
            _task_local_assembler(
                AssemblerType,
                buffer_range,
                assembler._matbuffer,
                assembler._rowbuffer,
                assembler._colbuffer,
                assembler._row_nalldofs,
                assembler._col_nalldofs,
            ),
        )
    end
    return assemblers
end

"""
    parallel_matrix_assembly(
        femms::AbstractArray{FEMM, 1},
        assemblers::AbstractVector{AT},
        matrixcomputation!::F,
        kind = :threaded,
    ) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, F<:Function}

Execute the assembly in parallel.
"""
function parallel_matrix_assembly(
    femms::AbstractArray{FEMM,1},
    assemblers::AbstractVector{AT},
    matrixcomputation!::F,
    kind = :threaded,
) where {FEMM<:AbstractFEMM,AT<:AbstractSysmatAssembler,F<:Function}
    length(femms) == length(assemblers) || error("Domains and assemblers mismatched")
    for j in eachindex(assemblers)
        assemblers[j]._buffer_pointer = 1
        setnomatrixresult(assemblers[j], true)
    end
    if kind == :threaded
        Threads.@threads for j in eachindex(femms)
            matrixcomputation!(femms[j], assemblers[j])
        end
    else
        Threads.@sync begin
            for j in eachindex(femms)
                Threads.@spawn let
                    matrixcomputation!(femms[j], assemblers[j])
                end
            end
        end
    end
    return true
end

end # module
