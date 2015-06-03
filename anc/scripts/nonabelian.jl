δE(1,11,1, linspace(-π/11, π/11, 20),linspace (-π, π, 20),0.02)

# -> Float64
# calculates the average (over specified points) non-abelian energy correction to n(th) band
# Note: MBZ:: (px, py) -> linspace(-π/q, π/q, N), linspace (-π, π, N)
δE(n::Int,q::Int,p::Int, px::Array{Float64, 1},py::Array{Float64, 1},κ::Float64) = mean([δE(n,q,p, x,y,κ) for y in py, x in px])

# -> Float64
# calculates non-abelian energy correction to n(th) band at position (kₓ⁰,ky) in the MBZ
δE(n::Int,q::Int,p::Int, k₀x::Float64,ky::Float64,κ::Float64) =  κ/2*sum([n′ != n ? norm(A(n,n′,q,p, k₀x,ky))^2 : 0.0 for n′ in 1:q])


# -> [Complex{Float64}, Complex{Float64}]
# calculates Berry connection A (vector quantity) at position (kₓ⁰,ky) in the MBZ for bands n and n′
function A(n::Int,n′::Int,q::Int,p::Int, k₀x::Float64,ky::Float64)
    # initializing hamiltonian matrices
    for M = (:H, :∇Hx, :∇Hy)
        @eval ($M) = zeros(Complex{Float64}, ($q,$q))
    end 

    # top right corner
    H[1,q] = -exp(-im*q*ky)
    ∇Hy[1,q] = im*q*exp(-im*q*ky) # ∂ky

    #bottom left corner
    H[q,1] = -exp(im*q*ky)
    ∇Hy[q,1] = -im*q*exp(im*q*ky) # ∂ky


    #upper subdiagonal
    ius = diagind(H,1)
    #lower subdiagonal
    ils = diagind(H,-1)
    
    for M = (:H, :∇Hx, :∇Hy), idx = (ius, ils)
        @eval ($M)[$idx] = -one(Complex{Float64})
    end 
    
    #main diagonal
    α = p/q
    for j in 1:q
        H[j,j] =  -2*cos(k₀x + 2*π*α*j)
      ∇Hx[j,j] =   2*sin(k₀x + 2*π*α*j)
    end 
        

    # diagonalize HH Hamiltonian -> E's and u's (eigs and eigvs)
    F = eigfact(H)
    U = F[:vectors]
    E = F[:values]

    # calculate denominator Float64
    denominator = E[n′] - E[n] 
    
    # calculate expectation value along kₓ (xnumerator)
    xnumerator = dot(U[:,n], ∇Hx*U[:,n′])
    # calculate expectation value along ky (ynumerator)
    ynumerator = dot(U[:,n], ∇Hy*U[:,n′])

    return [im * xnumerator / denominator, im * ynumerator / denominator]
end 



#import HH
using Base.Test


function foo(px::Array{Float64, 1}, py::Array{Float64, 1})
    Q = Array(Float64, length(py), length(px))

    for (j,x) in enumerate(px), (i,y) in enumerate(py)
        Q[i,j] = x-y
    end 
    Q
end 

bazz(px::Array{Float64, 1}, py::Array{Float64, 1}) = [x - y for y in py, x in px]


R = foo(linspace(1,10,10), linspace(1,5,5))
S = bazz(linspace(1,10,10), linspace(1,5,5))

@test_approx_eq R S
