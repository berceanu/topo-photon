module HH

import BP

##

#in-place builder of the q x q HH matrix in mom. sp.
function momsphhmat!(A::Matrix{Complex{Float64}},
                     kx0::Float64,ky::Float64,
                     p::Int)

    q::Int = size(A,1)
    
    A[1,q] = -exp(-im*q*ky)
    A[q,1] = -exp(im*q*ky)

    #upper subdiagonal
    for i in diagind(A,1)
        A[i] = -one(Complex{Float64})
    end 

    #lower subdiagonal
    for i in diagind(A,-1)
        A[i] = -one(Complex{Float64})
    end 

    #main diagonal
    α = p/q
    for (j,i) in enumerate(diagind(A))
        A[i] = -2*cos(kx0 + 2*π*α*j)
    end 

    nothing
end

momsphhmat!(A::Matrix{Complex{Float64}}, kx0::Float64,ky::Float64) = momsphhmat!(A, kx0,ky, 1)

function hhladder(p::Int, q::Int)
    ky = linspace(-π, π, 100)
    E = Array(Float64, (length(ky),q))
    H = zeros(Complex{Float64}, (q,q))
    for (i,k) in enumerate(ky)
        momsphhmat!(H, 0.,k, p)
        E[i,:] = eigvals(H)
    end
    return E
end

hhladder(q::Int) = hhladder(1,q)



function hhladder!(E::Matrix{Float64}, p::Int)

    ky = linspace(-π, π, size(E,1))
    q = size(E,2)
    H = zeros(Complex{Float64}, (q,q))
    for (i,k) in enumerate(ky)
        momsphhmat!(H, 0.,k, p)
        E[i,:] = eigvals(H)
    end
    nothing
end

hhladder!(E::Matrix{Float64}) = hhladder!(E, 1)

#finding energy difference
function ηzpe(N::Int, p::Int, q::Int, κ::Float64)
    α = p/q
    E = Array(Float64, (100,q))
    hhladder!(E, p)
    e0g = mean(E[:,1]) + .5*κ/(2π*α)

    nz = BP.countnonzeros(N)
    ft = ((n,m,a) -> one(Complex{Float64}), (n,m,a) -> one(Complex{Float64}),
          (n,m,a) -> exp(-im*2π*a*m), (n,m,a) -> exp(im*2π*a*m))
    ftex =  map(f -> (x, y, z) -> -f(x, y, z), ft)
    M = BP.genspmat(ftex...,(n,m,a) -> 1/2*κ*(n^2+m^2) + zero(Complex{Float64}), N,nz,α)

    e0o = real(eigs(M, nev=1, which=:SR, ritzvec=false)[1][1])

    return (e0o - e0g)/ (.5*κ/(2π*α))
end

ηzpe(p::Int, q::Int, κ::Float64) = ηzpe(35, p, q, κ)
ηzpe(q::Int, κ::Float64) = ηzpe(35, 1, q, κ)

#TODO: investigate Hermitian(), ishermitian().. may no longer need real()
#TODO: calculate only the first energy level of the q x q HH mat
#TODO: implement data type which gives Vector output for fixed κ OR q.

## 1. for each of the 100 ky values in [-π,­π], build the
## corresponding q x q HH matrix (which depends on k)
## and diagonalize it. we get q energy levels E(ky, 1..q) for each ky point.
## 2. take the mean (over all ky) of the *first* energy level E(ky,1)
## 3. generate the sparse matrix corresponding to HH + trap (no dissipation)
## 4. take its *first* eigenvalue
## 5. do the difference of the energies calculated at points 2. and 4.


##################################################
##################################################

# Computes the Harper-Hofstadter Hamiltonian matrix in momentum space
function momsphhmat(kx0::Float64, ky::Float64, α::Float64,q::Int)
    du = ones(Complex{Float64}, q-1) #upper diagonal
    d = Complex{Float64}[2*cos(kx0 + 2*π*α*j) for j in 1:q] #main diagonal
    mat = full(SymTridiagonal(d, du))
    mat[1,q] = exp(-im*q*ky)
    mat[q,1] = exp(im*q*ky)
    return -mat
end
# Computes the q energy levels E(p_y)
function hhladder(α::Float64, q::Int)
    kx0 = 0.
    ky = linspace(-π, π, 100)
    E = Array(Float64, 100,q)
    for c in 1:100
        M = Hermitian(momsphhmat(kx0, ky[c], α, q))
        E[c,:] = eigvals(M)
    end
    E
end

# level error
function ηl(q::Int,κ::Float64; N=35)
    α=1/q

    # generate Hamiltonian for HH+trap
    nonz = BP.countnonzeros(N)
    ft = ((n,m) -> -one(Complex{Float64}), (n,m) -> -one(Complex{Float64}), (n,m) -> -exp(-im*2π*α*m), (n,m) -> -exp(im*2π*α*m), (n,m) -> 1/2*κ*(n^2+m^2)+zero(Complex{Float64}))
    x = full(BP.genspmat(ft..., N,nonz))
    
    # get first 2 levels of spectrum of HH+trap
    spectrum = eigvals(Hermitian(x), 1:2)   

    #calculate level spacing
    2π*(α/κ) * (spectrum[2] - spectrum[1]) - 1
end


#finding bandwidth
#hhladder
function bwidth(q::Int)
    α=1/q
    E = hhladder(α,q)
    v = E[:,1]
    maximum(v) - minimum(v)
end
#finding bandgap
#hhladder
function bgap(q::Int)
    α=1/q
    E = hhladder(α,q)
    e1 = mean(E[:,1])
    e2 = mean(E[:,2])
    e2 - e1
end

end #module


#finding relative position
## function relpos(q::Int,κ::Float64; N=35)
##     α=1/q
##     E = hhladder(α,q)
##     e1 = mean(E[:,1])
##     e2 = mean(E[:,2])
##     b = e2 - e1
    
##     nonz = BP.countnonzeros(N)
##     ft = ((n,m) -> -one(Complex{Float64}), (n,m) -> -one(Complex{Float64}), (n,m) -> -exp(-im*2π*α*m), (n,m) -> -exp(im*2π*α*m), (n,m) -> 1/2*κ*(n^2+m^2)+zero(Complex{Float64}))
##     x = full(BP.genspmat(ft..., N,nonz))
##     e0o = eigvals(Hermitian(x), 1:1)[1]   
##     a = e0o - e1
    
##     a/b*100
## end
# Usage:
# map(relpos, qs, κò)