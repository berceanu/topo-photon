import BP

type WaveFunction
    ⅈ::Float64
    momφ::Matrix{Float64}
    spcφ::Matrix{Float64}
    momabs::Matrix{Float64}
    spcabs::Matrix{Float64}

    function WaveFunction(N,ft,ω,γ,κ,p,q,P,Np)
        α =  p/q
        kx = linspace(-π, π, Np)
        ky = linspace(-π, π, Np)
        nz = BP.countnonzeros(N)
        S = BP.genspmat(ft...,(n,m,a) -> ω + im*γ - 1/2*κ*(n^2+m^2), N,nz,α)
        X = reshape(S\P, N,N)
        Xk = BP.myfft2(X, kx, ky)
        new(sum(abs2(X)), abs(X), angle(X), abs(Xk), angle(Xk))
    end

end


ψ = WaveFunction(45,
                 ( 
                  (n,m,a) -> one(Complex{Float64}), (n,m,a) -> one(Complex{Float64}),
                  (n,m,a) -> exp(-im*2π*a*m), (n,m,a) -> exp(im*2π*a*m)
                  ),
                 -2., 0.001, 0.02, 1,11,
                 BP.δpmp(45; n0=5, m0=5),
                 200)





    



    
