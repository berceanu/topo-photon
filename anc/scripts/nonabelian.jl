import HH
import BP


# calculates zero-point energy, with and without nonabelian correction
ηzpenew(q::Int, κ::Float64) = (α=1/q; 4π*α/κ * endiffnonab(q, κ))
ηzpeold(q::Int, κ::Float64) = (α=1/q; 4π*α/κ * endiff(q, κ))


# calculates energy difference between numerical (exact) and theoretical energies, with the nonab corr.
endiffnonab(q::Int, κ::Float64) = er(q,κ) - (et(q,κ) + δE(q,κ))

# calculates energy difference between numerical (exact) and theoretical energies (without nonab corr)
endiff(q::Int, κ::Float64) =  er(q,κ) - et(q,κ)



et(q::Int,κ::Float64) =  et(q,1,κ)

# -> Float64
# calculates theoretical energy, without non-abelian correction
function et(q::Int,p::Int,κ::Float64)
    α = p/q

    gs = Array(Float64, 25, 25)
    HH.hhgrstate!(gs, p, q)

    e1 = mean(gs)

    e1 + 1/2*κ/(2π*α)
end 


er(q::Int,κ::Float64) = er(q,1,15,κ)

# -> Float64
# calculates numerical (exact) energy
function er(q::Int,p::Int,N::Int,κ::Float64)
    M = spzeros(Complex{Float64}, N^2,N^2)
    
    α = p/q
    BP.buildham_exact!(M, N,α,κ)

    real(eigs(M, nev=1, which=:SR, ritzvec=false)[1][1])
end


# calculates average over MBZ of the non-abelian correction to the 1(st) band for p=1 for certain trap
δE(q::Int,κ::Float64) = δE(1,q,1, linspace(-π/q, π/q, 20),linspace (-π, π, 20), κ)


# -> Float64
# calculates the average (over specified points) non-abelian energy correction to n(th) band
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
    H[1,q] = -exp(-im*ky)
    ∇Hy[1,q] = im*exp(-im*ky) # ∂ky

    #bottom left corner
    H[q,1] = -exp(im*ky)
    ∇Hy[q,1] = -im*exp(im*ky) # ∂ky


    #upper subdiagonal
    ius = diagind(H,1)
    #lower subdiagonal
    ils = diagind(H,-1)
    
    H[ius] = -exp(im*ky)
    H[ils] = -exp(-im*ky)

    ∇Hy[ius] = -im*exp(im*ky)
    ∇Hy[ils] = im*exp(-im*ky)

    
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



using Base.Test


## function foo(px::Array{Float64, 1}, py::Array{Float64, 1})
##     Q = Array(Float64, length(py), length(px))
##     for (j,x) in enumerate(px), (i,y) in enumerate(py)
##         Q[i,j] = x-y
##     end 
##     Q
## end 
## bazz(px::Array{Float64, 1}, py::Array{Float64, 1}) = [x - y for y in py, x in px]

## R = foo(linspace(1,10,10), linspace(1,5,5))
## S = bazz(linspace(1,10,10), linspace(1,5,5))

## @test_approx_eq R S

@test_approx_eq δE(1,5,1, 0.,0.,2.) 0.5378296443127227

v1 = A(1,2,5,1, 0.,0.);
v2 = Complex{Float64}[-0.0-0.518207im, -0.518207+0.0im];

@test_approx_eq_eps v1 v2 1e-6

#plotting

qs = 5:20;
y1 = [ηzpeold(q, 0.02)::Float64 for q in qs];
y2 = [ηzpenew(q, 0.02)::Float64 for q in qs];

y3 = HH.ηzpe(qs,0.02)
@test_approx_eq_eps y1 y3 1e-6

using PyPlot

# matplotlib parameters
matplotlib["rcParams"][:update](["axes.labelsize" => 22,
                                 "axes.titlesize" => 20,
                                 "font.size" => 18,
                                 "legend.fontsize" => 14,
                                 "axes.linewidth" => 1.5,
                                 "font.family" => "serif",
                                 "font.serif" => "Computer Modern Roman",
                                 "xtick.labelsize" => 20,
                                 "xtick.major.size" => 5.5,
                                 "xtick.major.width" => 1.5,
                                 "ytick.labelsize" => 20,
                                 "ytick.major.size" => 5.5,
                                 "ytick.major.width" => 1.5,
                                 "text.usetex" => true,
                                 "figure.autolayout" => true])


fig, ax = plt.subplots(figsize=(8, 3))
#TODO: use \text{} for plot subscripts. works only if nothing follows \text 

ax[:plot](qs, y1, "black", marker="o", ls="dashed", label=L"$E_{ex} - E_{th}$")
ax[:plot](qs, y2, "black", marker="o", label=L"$E_{ex} - (E_{th} + \delta E)$")

ax[:set_ylim](-1.1, 1.2)
ax[:set_xlim](qs[1], qs[end])

ax[:set_xlabel](L"$q$")
ax[:set_ylabel](L"$\eta_{zpe}$")

ax[:legend](loc="lower right")

fig[:savefig]("../../figures/nonabcorr.svg", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(fig)
