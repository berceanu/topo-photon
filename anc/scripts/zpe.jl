using Base.Test

import HH

#reload("HH")


#@test_approx_eq_eps(HH.ηzpe(11, 0.02), 0.9368071126414523, 1e-5)
#HH.ηzpe(11, 0.02)



N=45
nz = BP.countentries(N)
ft = ((n,m,a) -> one(Complex{Float64}), (n,m,a) -> one(Complex{Float64}),
      (n,m,a) -> exp(-im*2π*a*m),       (n,m,a) -> exp(im*2π*a*m))
ftex =  map(f -> (x, y, z) -> -f(x, y, z), ft)
#N²xN² matrix
M = BP.genspmat(ftex...,(n,m,a) -> 1/2*0.02*(n^2+m^2) + zero(Complex{Float64}), N,nz,1/11)

A = spzeros(Complex{Float64}, N^2,N^2)
BP.buildham_exact!(A, N,1/11, 0.02)
@test_approx_eq_eps(M, A, 1e-5)

@time for i=1:10^3; BP.buildham_exact!(A, N,1/11,0.02); end 


## for sym in {:landau,:symmetric,:exact}
##     @eval function $(symbol(string("buildham_", sym)))($(symbol(string("args_",sym))))
