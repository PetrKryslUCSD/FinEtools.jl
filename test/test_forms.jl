module mbilform_dot_1
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 2, 3, 4

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(psi)
    G = innerproduct(femm, geom, psi)
    # @show v' * G * v
    @test abs(v' * G * v - (W*L*t)) / (W*L*t) <= 1.0e-5
    G = bilform_dot(femm, geom, psi, DataCache(LinearAlgebra.I(1)))
    @test abs(v' * G * v - (W*L*t)) / (W*L*t) <= 1.0e-5
    true
end
test()
nothing
end

module mbilform_diffusion_1
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 2, 3, 4

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(psi)
    K = bilform_diffusion(femm, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    @test abs(v' * K * v - (0.0)) / (W*L*t) <= 1.0e-5
    true
end
test()
nothing
end

module mbilform_diffusion_2
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 2, 3, 4
    a, b, c, d = (-0.1, +0.3, +0.4, -0.5)

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(reshape([a + b*fens.xyz[j, 1] + c*fens.xyz[j, 2] + d*fens.xyz[j, 3]  for j in eachindex(fens)], count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(psi)
    K = bilform_diffusion(femm, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    @test abs(v' * K * v - (b^2 + c^2 + d^2) * (W*L*t)) / (W*L*t) <= 1.0e-5
    true
end
test()
nothing
end

module mbilform_diffusion_3
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 2, 3, 4
    a, b, c, d = (-0.1, +0.3, +0.4, -0.5)

    fens, fes = H20block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(reshape([a + b*fens.xyz[j, 1] + c*fens.xyz[j, 2] + d*fens.xyz[j, 3]  for j in eachindex(fens)], count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 3)))
    v = gathersysvec(psi)
    K = bilform_diffusion(femm, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    @test abs(v' * K * v - (b^2 + c^2 + d^2) * (W*L*t)) / (W*L*t) <= 1.0e-5
    true
end
test()
nothing
end

module mbilform_diffusion_4
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nw = 3, 4
    a, b, c = (-0.1, +0.3, +0.4)

    fens, fes = Q8block(L, W, nl, nw)
    geom = NodalField(fens.xyz)
    psi = NodalField(reshape([a + b*fens.xyz[j, 1] + c*fens.xyz[j, 2]  for j in eachindex(fens)], count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 3), t))
    v = gathersysvec(psi)
    K = bilform_diffusion(femm, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(2))))
    @test abs(v' * K * v - (b^2 + c^2) * (W*L*t)) / (W*L*t) <= 1.0e-5
    true
end
test()
nothing
end

module mbilform_diffusion_5
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nw = 3, 4
    a, b, c = (-0.1, +0.3, +0.4)

    fens, fes = Q8block(L, W, nl, nw)
    geom = NodalField(fens.xyz)
    psi = NodalField(reshape([a + b*fens.xyz[j, 1] + c*fens.xyz[j, 2]  for j in eachindex(fens)], count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 3), t))
    v = gathersysvec(psi)
    K = bilform_diffusion(femm, geom, psi, DataCache(1.0))
    @test abs(v' * K * v - (b^2 + c^2) * (W*L*t)) / (W*L*t) <= 1.0e-5
    true
end
test()
nothing
end

module mdistributedl1
using FinEtools
using Test
function test()
    W = 1.1
    L = 12.0
    t = 4.32
    nl, nt, nw = 2, 3, 4

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    fi = ForceIntensity([11.0])
    F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0)) / 667 <= 1.0e-5
    fi = ForceIntensity(11.0)
    F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0)) / 667 <= 1.0e-5

    psi = NodalField(fill(1.0 + 1.0im, count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    fi = ForceIntensity([11.0im])
    F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0im)) / 667 <= 1.0e-5
    fi = ForceIntensity(11.0im)
    F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0im)) / 667 <= 1.0e-5
    true
end
end
using .mdistributedl1
mdistributedl1.test()

module mbilform_convection_1
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 6.1
    L = 12.0
    t = 2.32
    nl, nw = 3, 4
    a, b, c = (-0.1, +0.3, +0.4)
    u_x, u_y = (3.1, -2.7)

    fens, fes = Q8block(L, W, nl, nw)
    geom = NodalField(fens.xyz)
    q = NodalField(reshape([a + b*fens.xyz[j, 1] + c*fens.xyz[j, 2]  for j in eachindex(fens)], count(fens), 1))
    numberdofs!(q)
    psi = deepcopy(q)
    psi.values .= 1.0

    u = NodalField(hcat(fill(u_x, count(fens), 1), fill(u_y, count(fens), 1)))
    numberdofs!(u)

    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 3), t))
    K = bilform_convection(femm, geom, u, q, DataCache(1.0))
    Q = gathersysvec(q)
    U = gathersysvec(u)
    Psi = gathersysvec(psi)
    # @show size(K)
    # @show (b * u_x + c * u_y) * (W*L*t)
    # @show Psi' * K * Q
    @test abs(Psi' * K * Q - (b * u_x + c * u_y) * (W*L*t)) / (W*L*t) <= 1.0e-5
    true
end
test()
nothing
end


module mbilform_div_grad_1
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 11.1
    L = 12.0
    t = 7.32
    nl, nt, nw = 2, 3, 4
    u_x, u_y, u_z = (3.1, -2.7, -0.77)
    mu = 0.00133

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)

    u = NodalField(hcat(fill(u_x, count(fens), 1), fill(u_y, count(fens), 1), fill(u_z, count(fens), 1)))
    numberdofs!(u)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(u)
    G = bilform_div_grad(femm, geom, u, DataCache(mu))
    # @show v' * G * v
    @test abs(v' * G * v - 0.0) / (W*L*t) <= 1.0e-5
    true
end
test()
nothing
end

module mbilform_div_grad_2
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 11.1
    L = 12.0
    t = 7.32
    nl, nt, nw = 2, 3, 4
    a, b, c, d = (-0.33, 2/3, -1.67, 2/7)
    mu = 0.13377

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)

    u = NodalField(
        hcat(
            reshape([a + b*fens.xyz[j, 1] + c*fens.xyz[j, 2] + d*fens.xyz[j, 3]  for j in eachindex(fens)], count(fens), 1),
            reshape([b + c*fens.xyz[j, 1] + d*fens.xyz[j, 2] + a*fens.xyz[j, 3]  for j in eachindex(fens)], count(fens), 1),
            reshape([c + d*fens.xyz[j, 1] + a*fens.xyz[j, 2] + b*fens.xyz[j, 3]  for j in eachindex(fens)], count(fens), 1)
            ))
    numberdofs!(u)
    gradu = [b c d; c d a; d a b]
    gradu_symm = (gradu + gradu') / 2
    int_true = 2 * mu * (W*L*t) * sum(gradu_symm.^2)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(u)
    G = bilform_div_grad(femm, geom, u, DataCache(mu))
    # @show v' * G * v
    @test abs(v' * G * v - int_true) / (int_true) <= 1.0e-5
    true
end
test()
nothing
end


module mbilform_lin_elastic_1
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 11.1
    L = 12.0
    t = 7.32
    nl, nt, nw = 2, 3, 4
    u_x, u_y, u_z = (3.1, -2.7, -0.77)
    mu = 0.00133

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)

    u = NodalField(hcat(fill(u_x, count(fens), 1), fill(u_y, count(fens), 1), fill(u_z, count(fens), 1)))
    numberdofs!(u)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(u)
    C = diagm([2*mu, 2*mu, 2*mu, mu, mu, mu])
    G = bilform_lin_elastic(femm, geom, u, DeforModelRed3D, DataCache(C))
    # @show v' * G * v
    @test abs(v' * G * v - 0.0) / (W*L*t) <= 1.0e-5
    true
end
test()
nothing
end

module mbilform_lin_elastic_2
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 11.1
    L = 12.0
    t = 7.32
    nl, nt, nw = 2, 3, 4
    a, b, c, d = (-0.33, 2/3, -1.67, 2/7)
    mu = 0.13377

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)

    u = NodalField(
        hcat(
            reshape([a + b*fens.xyz[j, 1] + c*fens.xyz[j, 2] + d*fens.xyz[j, 3]  for j in eachindex(fens)], count(fens), 1),
            reshape([b + c*fens.xyz[j, 1] + d*fens.xyz[j, 2] + a*fens.xyz[j, 3]  for j in eachindex(fens)], count(fens), 1),
            reshape([c + d*fens.xyz[j, 1] + a*fens.xyz[j, 2] + b*fens.xyz[j, 3]  for j in eachindex(fens)], count(fens), 1)
            ))
    numberdofs!(u)
    gradu = [b c d; c d a; d a b]
    gradu_symm = (gradu + gradu') / 2
    int_true = 2 * mu * (W*L*t) * sum(gradu_symm.^2)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(u)
    C = diagm([2*mu, 2*mu, 2*mu, mu, mu, mu])
    G = bilform_lin_elastic(femm, geom, u, DeforModelRed3D, DataCache(C))
    # v' * G * v / 2
    @test abs(v' * G * v  - int_true) / (int_true) <= 1.0e-5
    true
end
test()
nothing
end


module mbilform_div_grad_3
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 11.1
    L = 12.0
    t = 7.32
    nl, nt, nw = 12, 33, 24
    a, b, c, d = (-0.33, 2/3, -1.67, 2/7)
    mu = 0.13377

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)

    u = NodalField(
        hcat(
            reshape([a + b*fens.xyz[j, 1]^2 + c*fens.xyz[j, 2]^3 + d*fens.xyz[j, 3]  for j in eachindex(fens)], count(fens), 1),
            reshape([b + c*fens.xyz[j, 1] + d*fens.xyz[j, 2]^2 + a*fens.xyz[j, 3]^4  for j in eachindex(fens)], count(fens), 1),
            reshape([c + d*fens.xyz[j, 1]^3 + a*fens.xyz[j, 2] + b*fens.xyz[j, 3]^2  for j in eachindex(fens)], count(fens), 1)
            ))
    numberdofs!(u)


    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(u)

    G1 = bilform_div_grad(femm, geom, u, DataCache(mu))
    C = diagm([2*mu, 2*mu, 2*mu, mu, mu, mu])
    G2 = bilform_lin_elastic(femm, geom, u, DeforModelRed3D, DataCache(C))
    # @show v' * G1 * v
    # @show v' * G2 * v
    @test abs(v' * G1 * v - v' * G2 * v) / (v' * G1 * v) <= 1.0e-5
    true
end
test()
nothing
end



