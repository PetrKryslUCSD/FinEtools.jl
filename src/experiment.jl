module mmmmmmmmmmmm
using Compat.Test



function jac!(J::Array{Float64, 2}, X::Array{Float64, 2}, conn::C, gradNparams::Array{Float64, 2}) where {C}
    n = size(gradNparams, 1)
    @inbounds for j = 1:size(J, 2)
        @inbounds for i = 1:size(J, 1)
            Ja = 0.0
            @inbounds for k = 1:n
                Ja += X[conn[k], i] * gradNparams[k, j]
            end
            J[i, j] = Ja
        end
    end
    return J
end

function jact!(J::Array{Float64, 2}, X::Array{Float64, 2}, conn::C, gradNparams::Array{Float64, 2}) where {C}
    n = size(gradNparams, 1)
    @inbounds for j = 1:size(J, 2)
        @inbounds for i = 1:size(J, 1)
            Ja = 0.0
            @inbounds for k = 1:n
                Ja += X[i, conn[k]] * gradNparams[k, j]
            end
            J[i, j] = Ja
        end
    end
    return J
end

function jactup!(J::Array{Float64, 2}, X::Vector{NTuple{2,Float64}}, conn::C, gradNparams::Vector{NTuple{2,Float64}}) where {C}
    n = size(gradNparams, 1)
    for k = 1
        ck = conn[k]
        for j = 1:size(J, 2)
            @inbounds for i = 1:size(J, 1)
                J[i, j] = X[ck][i] * gradNparams[k][j]
            end
        end
    end
    for k = 2:n
        ck = conn[k]
        for j = 1:size(J, 2)
            @inbounds for i = 1:size(J, 1)
                J[i, j] += X[ck][i] * gradNparams[k][j]
            end
        end
    end
    return J
end


function ms!(J::Array{Float64, 2}, X::Array{Float64, 2}, Y::Array{Float64, 2}) 
    n = size(J, 1)
    @inbounds for j = 1:n
        @inbounds for i = 1:n
            Ja = 0.0
            @inbounds for k = 1:n
                Ja += X[k, i] * Y[k, j]
            end
            J[i, j] = Ja
        end
    end
    return J
end

function mst!(J::Array{Float64, 2}, X::Array{Float64, 2}, Y::Array{Float64, 2}) 
    n = size(J, 1)
    @inbounds for j = 1:n
        @inbounds for i = 1:n
            Ja = 0.0
            @inbounds for k = 1:n
                Ja += X[i, k] * Y[k, j]
            end
            J[i, j] = Ja
        end
    end
    return J
end

function test(N)
    gradNparams = rand(3,2)
    J = fill(0.0, 2, 2)
    X = rand(100000, 2)
    conn = (5, 7000, 19999)
    

    for i=1:N
        jac!(J, X, conn, gradNparams)
    end
end

function testt(N)
    gradNparams = rand(3,2)
    J = fill(0.0, 2, 2)
    X = rand(2, 100000)
    conn = (5, 7000, 19999)
    

    for i=1:N
        jact!(J, X, conn, gradNparams)
    end
end

function testtup(N)
    gradNparams = Vector{NTuple{2,Float64}}(3)
    J = fill(0.0, 2, 2)
    X = Vector{NTuple{2,Float64}}(100000)
    conn = (5, 7000, 19999)
    

    for i=1:N
        jactup!(J, X, conn, gradNparams)
    end
end

function testms(N)
    M = 500
    Y = rand(M,M)
    J = fill(0.0, M, M)
    X = rand(M, M)
    
    for i=1:N
        ms!(J, X, Y)
    end
end

function testmst(N)
    M = 500
    Y = rand(M,M)
    J = fill(0.0, M, M)
    X = rand(M, M)
    
    for i=1:N
        mst!(J, X, Y)
    end
end


function testmsb(N)
    M = 500
    Y = rand(M,M)
    J = fill(0.0, M, M)
    X = rand(M, M)
    
    for i=1:N
        A_mul_B!(J, X, Y)
    end
end

function testmstb(N)
    M = 500
    Y = rand(M,M)
    J = fill(0.0, M, M)
    X = rand(M, M)
    
    for i=1:N
        At_mul_B!(J, X, Y)
    end
end

function mv_product!(r::Vector{Float64}, Ke::Array{Float64, 2}, v::Vector{Float64})
    @assert size(Ke, 1) == size(Ke, 2)
    @assert length(v) == size(Ke, 2)
    Kedim = size(Ke, 1)
    mx::Int = 0
    nx::Int = 0
    @inbounds for mx = 1:Kedim
        accumulator::Float64 = 0.0
        @inbounds for nx = 1:Kedim
            accumulator += Ke[mx, nx] * v[nx]
        end
        r[mx] = accumulator
    end
    return true
end

function testmvp(N)
    M = 5000
    Y = rand(M)
    J = fill(0.0, M)
    X = rand(M, M)
    
    for i=1:N
        mv_product!(J, X, Y)
    end
end

function testmvpb(N)
    M = 5000
    Y = rand(M)
    J = fill(0.0, M)
    X = rand(M, M)
    
    for i=1:N
        A_mul_B!(J, X, Y)
    end
end

end
using .mmmmmmmmmmmm
# N = 100000000
# @time mmmmmmmmmmmm.test(N)
# @time mmmmmmmmmmmm.testt(N)
# @time mmmmmmmmmmmm.testtup(N)

# N = 10
# @time mmmmmmmmmmmm.testms(N)
# @time mmmmmmmmmmmm.testmst(N)
# @time mmmmmmmmmmmm.testmsb(N)
# @time mmmmmmmmmmmm.testmstb(N)

N = 100
@time mmmmmmmmmmmm.testmvpb(N)
@time mmmmmmmmmmmm.testmvp(N)