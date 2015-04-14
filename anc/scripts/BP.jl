module BP

using Polynomials


type WaveFunction
    N::Int
    int::Float64
    ψ::Vector{Complex{Float64}}
end

function WaveFunction(ω::Float64,P::Vector{Complex{Float64}},gauge::Symbol,
                      α::Float64,γ::Float64,κ::Float64)

    gauges = [:landau => 
              ((n,m,a) -> one(Complex{Float64}), (n,m,a) -> one(Complex{Float64}),
               (n,m,a) -> exp(-im*2π*a*m), (n,m,a) -> exp(im*2π*a*m)),
              :symmetric =>
              ((n,m,a) -> exp(-im*π*a*n), (n,m,a) -> exp(im*π*a*n),
               (n,m,a) -> exp(-im*π*a*m), (n,m,a) -> exp(im*π*a*m))]

    N::Int = sqrt(length(P))
    nz::Int = countnonzeros(N) #

    S = genspmat(gauges[gauge]...,(n,m,a) -> ω + im*γ - 1/2*κ*(n^2+m^2), N,nz,α)
    X = S\P

    return WaveFunction(N,sum(abs2(X)),X)
end

WaveFunction(ω::Float64,P::Vector{Complex{Float64}},gauge::Symbol) =
    WaveFunction(ω,P,gauge,1/11,0.001,0.1)

WaveFunction(ω::Float64,P::Vector{Complex{Float64}}) = WaveFunction(ω,P,:landau)

type Spectrum
    N::Int
    gauge::Symbol
    pump::Vector{Complex{Float64}}
    νs::Vector{Float64}
    intensity::Vector{Float64}
    states::Vector{WaveFunction}
end

function Spectrum(ν::Vector{Float64},P::Vector{Complex{Float64}},gauge::Symbol)
    statevec = Array(WaveFunction, length(ν))
    intvec = Array(Float64, length(ν))
    for (i,ω) in enumerate(ν)
        statevec[i] = WaveFunction(ω, P, gauge)
        intvec[i] = statevec[i].int
    end 

    N::Int = sqrt(length(P))
        
    return Spectrum(N, gauge, P, ν, intvec, statevec)
end

Spectrum(ν::Vector{Float64},P::Vector{Complex{Float64}}) = Spectrum(ν,P,:landau)

function getstate(s::Spectrum, ω::Float64)
    i::Int = indmin(abs(s.νs .- ω))
    return s.states[i]
end 
##

getm(i::Int64,N::Int64) = div(i-1,N)-div(N-1,2)
getn(i::Int64,N::Int64) = div(N-1,2)-rem(i-1,N)
geti(m::Int,n::Int,N::Int)=(m+div(N-1,2))*N+(div(N-1,2)-n)+1

countentries(N::Int) = N^2 + 8 + 4*(N-2)*3 + (N-2)^2*4

macro hambody(fself, fleft, fright, fup, fdown)
    return quote
        border::Int = div(N-1,2)
        for m in -border:border, n in -border:border
            i  = geti(m,n,N)
            S[i,i] = $fself
        end 
        for m in -border+1:border, n in -border:border
            i  = geti(m,n,N)
            S[i,i-N] = $fleft
        end 
        for m in -border:border-1, n in -border:border
            i  = geti(m,n,N)
            S[i,i+N] = $fright
        end 
        for m in -border:border, n in -border:border-1
            i  = geti(m,n,N)
            S[i,i-1] = $fup
        end 
        for m in -border:border, n in -border+1:border
            i  = geti(m,n,N)
            S[i,i+1] = $fdown
        end
    end 
end 

function buildham_landau!(S::SparseMatrixCSC{Complex{Float64},Int}, N::Int,α::Float64,κ::Float64,γ::Float64,ω::Float64)
    @hambody(ω + im*γ - 1/2*κ*(n^2+m^2), 1, 1, exp(-im*2π*α*m), exp(im*2π*α*m))
end

function buildham_exact!(S::SparseMatrixCSC{Complex{Float64},Int}, N::Int,α::Float64,κ::Float64)
    @hambody(1/2*κ*(n^2+m^2), -1, -1, -exp(-im*2π*α*m), -exp(im*2π*α*m))
end

function buildham_symmetric!(S::SparseMatrixCSC{Complex{Float64},Int}, N::Int,α::Float64,κ::Float64,γ::Float64,ω::Float64)
    @hambody(ω + im*γ - 1/2*κ*(n^2+m^2), exp(-im*π*a*n), exp(im*π*a*n), exp(-im*π*a*m), exp(im*π*a*m))
end



function genspmat(l::Function,r::Function,u::Function,d::Function,s::Function, N::Int,nz::Int,α::Float64)
    iseven(N) && throw(ArgumentError("invalid system size N=$N. N must be odd"))
    # Preallocate
    I = Array(Int64,nz)
    J = Array(Int64,nz)
    V = Array(Complex{Float64},nz)

    function setnzelem(i::Int,n::Int,m::Int; pos::ASCIIString = "self")
        if pos=="left"
            k += 1
            J[k] = i-N; I[k] = i; V[k] = l(n,m,α)
        elseif pos=="right"
            k += 1
            J[k] = i+N; I[k] = i; V[k] = r(n,m,α)
        elseif pos=="up"
            k += 1
            J[k] = i-1; I[k] = i; V[k] = u(n,m,α)
        elseif pos=="down"
            k += 1
            J[k] = i+1; I[k] = i; V[k] = d(n,m,α)
        elseif pos=="self"
            k += 1
            J[k] = i; I[k] = i; V[k] = s(n,m,α)
        end
    end
            
    # maximum value of m or n indices
    maxm = div(N-1,2)

    k = 0
    for i in 1:N^2
        m = getm(i,N)
        n = getn(i,N)
        setnzelem(i,n,m; pos="self")
        #corners
        #top left
        if n==maxm && m==-maxm
            setnzelem(i,n,m; pos="right")
            setnzelem(i,n,m; pos="down")
        #top right
        elseif n==maxm && m==maxm
            setnzelem(i,n,m; pos="left")
            setnzelem(i,n,m; pos="down")
        #bottom right
        elseif n==-maxm && m==maxm 
            setnzelem(i,n,m; pos="left")
            setnzelem(i,n,m; pos="up")
        #bottom left
        elseif n==-maxm && m==-maxm 
            setnzelem(i,n,m; pos="right")
            setnzelem(i,n,m; pos="up")
        #edges
        #top
        elseif n == maxm
            setnzelem(i,n,m; pos="right")
            setnzelem(i,n,m; pos="left")
            setnzelem(i,n,m; pos="down")
        #right
        elseif m == maxm
            setnzelem(i,n,m; pos="left")
            setnzelem(i,n,m; pos="up")
            setnzelem(i,n,m; pos="down")
        #bottom
        elseif n == -maxm
            setnzelem(i,n,m; pos="left")
            setnzelem(i,n,m; pos="up")
            setnzelem(i,n,m; pos="right")
        #left
        elseif m == -maxm
            setnzelem(i,n,m; pos="down")
            setnzelem(i,n,m; pos="up")
            setnzelem(i,n,m; pos="right")
        else #bulk
            setnzelem(i,n,m; pos="down")
            setnzelem(i,n,m; pos="up")
            setnzelem(i,n,m; pos="right")
            setnzelem(i,n,m; pos="left")
        end
    end

    return sparse(I,J,V)
end



########################
#various pumping schemes
########################
function δpmp(N::Int; A=1., seed=0, σ=0., n0=0, m0=0)
    i = (m0+div(N-1,2)) * N + (div(N-1,2)-n0) + 1
    f = zeros(Complex{Float64}, N^2)
    f[i] = A * one(Complex{Float64})
    f
end
function gausspmp(N::Int; A=1., seed=0, σ=1., n0=0, m0=0)
    x0=m0 
    y0=n0
    f = zeros(Complex{Float64}, N,N)
    m = [-div(N-1,2):div(N-1,2)]
    n = [div(N-1,2):-1:-div(N-1,2)]
    for c in 1:N, l in 1:N
        x = m[c]
        y = n[l]
        f[l,c] = A*exp(-1/(2σ^2)*((x-x0)^2 + (y-y0)^2))
    end
    reshape(f, N^2)
end
function randpmp(N::Int; A=1., seed=123, σ=0., n0=0, m0=0)
    # seed the RNG #
    srand(seed)
    # generate matrix of random phases in interval [0,2π)
    ϕ = 2π .* rand(N^2)
    A .* exp(im .* ϕ)
end
function homopmp(N::Int; A=1., seed=0, σ=0., n0=0, m0=0)
    A .* ones(Complex{Float64}, N^2)
end


#########################
#arbitrary resolution fft
#########################
function myfft2(ψr::Matrix{Complex{Float64}}, k1::Float64, k2::Float64, xs1::Float64, xs2::Float64, Δx1::Float64, Δx2::Float64)
    (N1,N2) = size(ψr)

    s = zero(Complex{Float64})
    for n2 in 1:N2, n1 in 1:N1
        xn1 = xs1 + (n2-1)*Δx1 #x
        xn2 = xs2 + (n1-1)*Δx2 #y
        cexp = exp(-im*(k1*xn1 + k2*xn2))
        s += ψr[n1,n2]*cexp
    end
    s
end
function myfft2(ψr::Matrix{Complex{Float64}}, k1, k2)
    N1 = length(k2); N2 = length(k1)
    

    out = Array(Complex{Float64}, N1,N2)
    for j in 1:N2, i in 1:N1
        out[i,j] = myfft2(ψr, k1[j], k2[i], 0., 0., 1., 1.)
    end
    out
end

# Magnetic Brillouin Zone
function mbz(data, q, N)
    #l = div(N-1, q)
    l = N

    V = zeros(Float64,size(data)[1],l)
    
    for i in 1:l
        idx = i
        while idx <= q*N+1
               V[:,i] += data[:,idx]
            idx += l
        end
    end

    V/(4π^2)
end

########################
#comparison to analytics
########################
function compute_hermite_polynomial(n)
    P = Poly([1])
    const x = Poly([0; 1])                                                                                 
    for i = 1:n
        P = 2x*P - polyder(P)
    end
    P
end
function χ(kx0, ky, α, β)
    l = sqrt(2π*α)

    sum = zero(Complex{Float64})
    for j in -20:20 #truncate the sum
        H = polyval(compute_hermite_polynomial(β), kx0/l + j*l)
        sum += exp(-im*ky*j) * exp(-(kx0 + j*l^2)^2/(2l^2)) * H
    end

    sum
end

end
