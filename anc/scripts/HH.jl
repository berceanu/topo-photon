module HH


import BP


#in-place builder of the q x q HH matrix in mom. sp.
momsphhmat!(A::Matrix{Complex{Float64}}, kx0::Float64,ky::Float64) = momsphhmat!(A, kx0,ky, 1)

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


hhladder(q::Int) = hhladder(1,q)
hhladder!(E::Matrix{Float64}) = hhladder!(E, 1)

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


# ground state energy of the HH hamiltonian
hhgrstate!(ve::Vector{Float64}, q::Int) = hhgrstate!(ve, 1, q)

function hhgrstate!(ve::Vector{Float64}, p::Int, q::Int)
    ky = linspace(-π, π, length(ve))

    H = zeros(Complex{Float64}, (q,q))
    for (i,k) in enumerate(ky)
        momsphhmat!(H, 0.,k, p)
        ve[i] = eigmin(H)
    end
    nothing
end


#finding energy difference
ηzpe(q::Int, κ::Float64) = ηzpe(1, q, κ)
ηzpe(p::Int, q::Int, κ::Float64) = (N=15; A = spzeros(Complex{Float64}, N^2,N^2); ηzpe(A, p, q, κ))


function ηzpe(M::SparseMatrixCSC{Complex{Float64},Int}, p::Int, q::Int, κ::Float64)
    N::Int = sqrt(size(M,1))
    α::Float64 = p/q
    gs = Array(Float64, 25)
    hhgrstate!(gs, p, q)
    e1 = mean(gs)
    et = e1 + 1/2*κ/(2π*α)

    BP.buildham_exact!(M, N,α,κ)

    er = real(eigs(M, nev=1, which=:SR, ritzvec=false)[1][1])

    return 4π*α/κ * (er - et)
end


ηzpe(q::Int, κs::Vector{Float64}) = vec(ηzpe([q], κs))

ηzpe(qs::UnitRange{Int}, κ::Float64) = ηzpe([qs], κ)
ηzpe(qs::Vector{Int}, κ::Float64) = vec(ηzpe(qs, [κ]))

ηzpe(qs::UnitRange{Int}, κs::Vector{Float64}) = ηzpe([qs], κs)
function ηzpe(qs::Vector{Int}, κs::Vector{Float64})
    N=15
    A = spzeros(Complex{Float64}, N^2,N^2)
    η = Array(Float64, length(qs), length(κs))
    for col=1:length(κs), row=1:length(qs)
        η[row,col] = ηzpe(A, 1, qs[row], κs[col])
    end
    return η
end 



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
    nonz = BP.countentries(N)
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
    
##     nonz = BP.countentries(N)
##     ft = ((n,m) -> -one(Complex{Float64}), (n,m) -> -one(Complex{Float64}), (n,m) -> -exp(-im*2π*α*m), (n,m) -> -exp(im*2π*α*m), (n,m) -> 1/2*κ*(n^2+m^2)+zero(Complex{Float64}))
##     x = full(BP.genspmat(ft..., N,nonz))
##     e0o = eigvals(Hermitian(x), 1:1)[1]   
##     a = e0o - e1
    
##     a/b*100
## end
# Usage:
# map(relpos, qs, κò)

## bw = [bwidth(q) for q in qs]
## bg = [bgap(q) for q in qs]

## hheig = HH.hhladder(α,q)
## ladder = [mean(hheig[:,1]) + (β + 1/2)*κ/(2π*α) for β in 0:14]
## ## using Color
## ## dc = distinguishable_colors(q)
## a = Layer[]
## for i = 1:q
##     push!(a, layer(x=linspace(-π, π, 100), y=hheig[:,i], Geom.line, Theme(default_color=dc[i]))[1] )
## end
## plot(a)
