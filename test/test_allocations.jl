module mbilform_dot_1
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 1.1
    L = 12.0
    t = 0.3
    ref = 20
    nl, nt, nw = ref*2, ref*3, ref*4

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(psi)
    @time G = innerproduct(femm, geom, psi)
    # @show v' * G * v
    @test abs(v' * G * v - (W*L*t)) / (W*L*t) <= 1.0e-5
    @time G = bilform_dot(femm, geom, psi, DataCache(LinearAlgebra.I(1)))
    @test abs(v' * G * v - (W*L*t)) / (W*L*t) <= 1.0e-5


    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(psi)
    @time K = bilform_diffusion(femm, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
   
    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    a, b, c, d = (-0.1, +0.3, +0.4, -0.5)
    psi = NodalField(reshape([a + b*fens.xyz[j, 1] + c*fens.xyz[j, 2] + d*fens.xyz[j, 3]  for j in eachindex(fens)], count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(psi)
    @time K = bilform_diffusion(femm, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 3)))
    v = gathersysvec(psi)
    @time K = bilform_diffusion(femm, geom, psi, DataCache(1.0))
    true
end
test()
function test2()
    W = 1.1
    L = 12.0
    t = 4.32
    ref = 20
    nl, nt, nw = ref*2, ref*3, ref*4

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    fi = ForceIntensity([11.0])
    @time F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0)) / 667 <= 1.0e-5
    fi = ForceIntensity(11.0)
    @time F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0)) / 667 <= 1.0e-5

    psi = NodalField(fill(1.0 + 1.0im, count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    fi = ForceIntensity([11.0im])
    @time F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0im)) / 667 <= 1.0e-5
    fi = ForceIntensity(11.0im)
   @time  F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0im)) / 667 <= 1.0e-5
    true
end
test2()
function test3()
    W = 6.1
    L = 12.0
    t = 2.32
    ref = 20
    nl, nt, nw = ref*2, ref*3, ref*4
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
    @time K = bilform_convection(femm, geom, u, q, DataCache(1.0))
    Q = gathersysvec(q)
    U = gathersysvec(u)
    Psi = gathersysvec(psi)
    # @show size(K)
    # @show (b * u_x + c * u_y) * (W*L*t)
    # @show Psi' * K * Q
    @test abs(Psi' * K * Q - (b * u_x + c * u_y) * (W*L*t)) / (W*L*t) <= 1.0e-5
    true
end
test3()
function test4()
    W = 11.1
    L = 12.0
    t = 7.32
    ref = 20
    nl, nt, nw = ref*2, ref*3, ref*4
    u_x, u_y, u_z = (3.1, -2.7, -0.77)
    mu = 0.00133

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)

    u = NodalField(hcat(fill(u_x, count(fens), 1), fill(u_y, count(fens), 1), fill(u_z, count(fens), 1)))
    numberdofs!(u)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(u)
    @time G = bilform_div_grad(femm, geom, u, DataCache(mu))
    # @show v' * G * v
    @test abs(v' * G * v - 0.0) / (W*L*t) <= 1.0e-5
    true
end
test4()
function test5()
    W = 11.1
    L = 12.0
    t = 7.32
    ref = 20
    nl, nt, nw = ref*2, ref*3, ref*4
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
    @time G = bilform_div_grad(femm, geom, u, DataCache(mu))
    # @show v' * G * v
    
    true
end
test5()
function test6()
    W = 11.1
    L = 12.0
    t = 7.32
    ref = 20
    nl, nt, nw = ref*2, ref*3, ref*4
    u_x, u_y, u_z = (3.1, -2.7, -0.77)
    mu = 0.00133

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)

    u = NodalField(hcat(fill(u_x, count(fens), 1), fill(u_y, count(fens), 1), fill(u_z, count(fens), 1)))
    numberdofs!(u)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(u)
    C = diagm([2*mu, 2*mu, 2*mu, mu, mu, mu])
    @time G = bilform_lin_elastic(femm, geom, u, DeforModelRed3D, DataCache(C))
    # @show v' * G * v
    
    true
end
test6()
function test7()
    W = 11.1
    L = 12.0
    t = 7.32
    ref = 20
    nl, nt, nw = ref*2, ref*3, ref*4
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
    @time G = bilform_lin_elastic(femm, geom, u, DeforModelRed3D, DataCache(C))
    # v' * G * v / 2
    
    true
end
test7()
function test8()
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

    @time G1 = bilform_div_grad(femm, geom, u, DataCache(mu))
    C = diagm([2*mu, 2*mu, 2*mu, mu, mu, mu])
    @time G2 = bilform_lin_elastic(femm, geom, u, DeforModelRed3D, DataCache(C))
    
    true
end
test8()

function test9()
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

    @time G1 = bilform_div_grad(femm, geom, u, DataCache(mu))
    # @test G1 - G1' == spzeros(size(G1)...)
    C = diagm([2*mu, 2*mu, 2*mu, mu, mu, mu])
    @time G2 = bilform_lin_elastic(femm, geom, u, DeforModelRed3D, DataCache(C))
    
    true
end
test9()

function test10()
    W = 11.1
    L = 12.0
    t = 7.32
    nl, nt, nw = 12, 33, 24
    mu = 0.13377

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)

    u = NodalField(
        hcat(
            reshape([ fens.xyz[j, 1] + 3 * fens.xyz[j, 2] for j in eachindex(fens)], count(fens), 1),
            reshape([ fens.xyz[j, 2] - 2 * fens.xyz[j, 3] for j in eachindex(fens)], count(fens), 1),
            reshape([ fens.xyz[j, 1] - 2 * fens.xyz[j, 3] for j in eachindex(fens)], count(fens), 1)
            ))

    numberdofs!(u)


    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(u)

    @time G1 = bilform_div_grad(femm, geom, u, DataCache(mu))
    # @test G1 - G1' == spzeros(size(G1)...)
    C = diagm([2*mu, 2*mu, 2*mu, mu, mu, mu])
    @time G2 = bilform_lin_elastic(femm, geom, u, DeforModelRed3D, DataCache(C))
    
    true
end
test10()

function test11()
    W = 11.1
    L = 12.0
    t = 7.32
    nl, nt, nw = 12, 33, 24
    mu = 0.13377

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)

    u = NodalField(
            hcat(
                reshape([ fens.xyz[j, 2] for j in eachindex(fens)], count(fens), 1),
                reshape([-fens.xyz[j, 1] for j in eachindex(fens)], count(fens), 1),
                reshape([ fens.xyz[j, 2]^2 for j in eachindex(fens)], count(fens), 1)
                ))

    numberdofs!(u)


    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(u)

    @time G1 = bilform_div_grad(femm, geom, u, DataCache(mu))
    # @test G1 - G1' == spzeros(size(G1)...)
    C = diagm([2*mu, 2*mu, 2*mu, mu, mu, mu])
    @time G2 = bilform_lin_elastic(femm, geom, u, DeforModelRed3D, DataCache(C))
    
    true
end
test11()
nothing
end
