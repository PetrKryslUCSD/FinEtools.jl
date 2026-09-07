module testerror1
using FinEtools
using Test

function test()
    fens, fes = T3block(1.0, 1.0, 1, 1)
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 1))
    exact(x) = 1.0
    integdomain = IntegDomain(fes, TriRule(3))
    femm = FEMMBase(integdomain)
    err = L2error(femm, geom, u, exact)
    @test size(err.values) == (2, 1)
    @test err.values[1] +  err.values[2] ≈ sqrt(2) atol = 1e-8
    return true
end

end

using .testerror1
testerror1.test()

module testerror2
using FinEtools
using Test

function test()
    fens, fes = T3block(1.0, 1.0, 1, 1)
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2))
    exact(x) = [1.0,1.0]
    integdomain = IntegDomain(fes, TriRule(3))
    femm = FEMMBase(integdomain)
    err = L2error(femm, geom, u, exact)
    @test size(err.values) == (2, 1)
    @test err.values[1] +  err.values[2] ≈ 2 atol = 1e-8
    return true
end

end

using .testerror2
testerror2.test()

module testerror3
using FinEtools
using Test

function test()
    fens, fes = T4block(1.0, 1.0, 1.0, 1, 1, 1)
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3))
    exact(x) = [1.0,1.0,1.0]
    integdomain = IntegDomain(fes, TetRule(4))
    femm = FEMMBase(integdomain)
    err = L2error(femm, geom, u, exact)
    @test size(err.values) == (6, 1)
    @test sum(err.values) ≈ 3*sqrt(2) atol = 1e-8
    return true
end

end

using .testerror3
testerror3.test()