module mparallelassembly_1
using FinEtools
using FinEtools.ParallelAssemblyModule: make_assembler, make_task_assemblers, parallel_matrix_assembly, make_matrix!
using LinearAlgebra
using Test

function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 12, 33, 24
    nl, nt, nw = 2, 3, 2
    ntasks = 2

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    nl = collect(1:3)
    setebc!(psi, nl, true, ones(Int, length(nl)), 0.0)
    numberdofs!(psi)
    v_f = gathersysvec(psi)

    chunks = Iterators.partition(eachindex(fes), Int(round(count(fes) / ntasks)))

    femms = FEMMBase[]
    for ch in chunks
        push!(femms, FEMMBase(IntegDomain(subset(fes, ch), GaussRule(3, 2))))
    end

    function matrixcomputation!(femm, assembler)
        bilform_diffusion(femm, assembler, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    end

    assembler = make_assembler(femms, SysmatAssemblerSparse, psi)
    assemblers = make_task_assemblers(femms, assembler, SysmatAssemblerSparse, psi)
    parallel_matrix_assembly(femms, assemblers, matrixcomputation!)
    K = make_matrix!(assembler)

    K_ff = matrix_blocked_ff(K, nfreedofs(psi))
    result = abs(v_f' * K_ff * v_f)

    ass = SysmatAssemblerFFBlock(nfreedofs(psi))
    K_ff2 = bilform_diffusion(FEMMBase(IntegDomain(fes, GaussRule(3, 2))), ass, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    @test norm(K_ff - K_ff2) / norm(K_ff2) <= 1.0e-5
    @test abs(v_f' * K_ff2 * v_f - (result)) / (result) <= 1.0e-5
    true
end
test()
nothing
end

module mparallelassembly_2
using FinEtools
using FinEtools.ParallelAssemblyModule: make_assembler, make_task_assemblers, parallel_matrix_assembly, make_matrix!
using LinearAlgebra
using Test

function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 12, 33, 24
    nl, nt, nw = 2, 3, 2
    ntasks = 2

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    nl = collect(1:3)
    setebc!(psi, nl, true, ones(Int, length(nl)), 0.0)
    numberdofs!(psi)
    v_f = gathersysvec(psi)

    chunks = Iterators.partition(eachindex(fes), Int(round(count(fes) / ntasks)))

    femms = FEMMBase[]
    for ch in chunks
        push!(femms, FEMMBase(IntegDomain(subset(fes, ch), GaussRule(3, 2))))
    end

    function matrixcomputation!(femm, assembler)
        bilform_diffusion(femm, assembler, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    end

    assembler = make_assembler(femms, SysmatAssemblerSparseSymm, psi)
    start_assembler!(assembler)
    assemblers = make_task_assemblers(femms, assembler, SysmatAssemblerSparseSymm, psi)
    parallel_matrix_assembly(femms, assemblers, matrixcomputation!)
    K = make_matrix!(assembler)

    K_ff = matrix_blocked_ff(K, nfreedofs(psi))
    result = abs(v_f' * K_ff * v_f)

    ass = SysmatAssemblerFFBlock(nfreedofs(psi))
    K_ff2 = bilform_diffusion(FEMMBase(IntegDomain(fes, GaussRule(3, 2))), ass, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    @test norm(K_ff - K_ff2) / norm(K_ff2) <= 1.0e-5
    @test abs(v_f' * K_ff2 * v_f - (result)) / (result) <= 1.0e-5
    true
end
test()
nothing
end


