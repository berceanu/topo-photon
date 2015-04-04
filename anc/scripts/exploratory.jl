import BP
using DSP

##################################################
#system parameters
sp=1
sq=11
sγ=0.001
sκ=0.02
sN=35
n=5
m=5
tbσ=5.0
pmptype=BP.δpmp
#pmptype=BP.gausspmp
#pmptype=BP.homopmp
#pmptype=BP.randpmp
sω1=-3.5
sω2=-2.4
dδ=0.001
ft = (
    (n,m,a) -> one(Complex{Float64}), (n,m,a) -> one(Complex{Float64}),
    (n,m,a) -> exp(-im*2π*a*m), (n,m,a) -> exp(im*2π*a*m))
#Symmetric gauge
## ft = (
##     (n,m,a) -> exp(-im*π*a*n), (n,m,a) -> exp(im*π*a*n),
##     (n,m,a) -> exp(-im*π*a*m), (n,m,a) -> exp(im*π*a*m))
sω0=-2.96101
sζ=10000
##################################################

ftex =  map(f -> (x, y, z) -> -f(x, y, z), ft);

g =  gcd(sp,sq)
p =  div(sp,g)
q =  div(sq,g)
α =  p/q

P = pmptype(sN; A=1., seed=1234, σ=tbσ, n0=n, m0=m)

freq =  sω1:dδ:sω2
println("datapoints in spectrum: $(length(freq))")

nz = BP.countnonzeros(sN)

function spect(ν::FloatRange{Float64})
    l = length(ν)
    v = Array(Float64, l)
    for i in 1:l
        S = BP.genspmat(ft...,(n,m,a) -> ν[i] + im*sγ - 1/2*sκ*(n^2+m^2), sN,nz,α)
        X = S\P
        v[i] =  sum(abs2(X))
    end
    v
end

spect(0.:0.)
ⅈ = spect(freq)

# we need to convert from sparse to dense in order to use eigvals
M = full(BP.genspmat(ftex...,(n,m,a) -> 1/2*sκ*(n^2+m^2) + zero(Complex{Float64}),sN,nz,α))

#M is the "exact" hamiltonian matrix, without dissipation and pumping
#exact spectrum, first 15 eigenvalues
spectrum =  eigvals(Hermitian(M), 1:15);

η =  indmin(abs(sω0-spectrum))
β =  η-1

Np =  q*sN+1
ky =  2π*DSP.fftshift(DSP.fftfreq(sN))
kx =  2π*DSP.fftshift(DSP.fftfreq(Np));

#returns the nth eigenvector
F = reshape(eigfact(Hermitian(M), η:η)[:vectors], sN,sN)
#exact wf
ψrex = abs2(F);
ψkex = abs2(BP.myfft2(F, kx, ky))
ψkmbzex = BP.mbz(ψkex, q,sN);


#matrix of linear system with dissipation
S = BP.genspmat(ft..., (n,m,a) -> sω0 + im*sγ - 1/2*sκ*(n^2+m^2), sN,nz,α)
X =   reshape(S\P, sN,sN);
ψr =   abs2(X)
ψk =   abs2(BP.myfft2(X, kx, ky))
ψkmbz =   BP.mbz(ψk, q,sN);

function getχ()
    yaxis = linspace(-π,π,sN)
    v = Array(Complex{Float64}, sN)
    for i = 1:sN
        v[i] = BP.χ(0.,yaxis[i],α,β)
    end
    radical = sqrt(sqrt(2/q) / (2π*2^β * factorial(β) * 2π*α))
    radical .* v
end

χ = getχ()

xm =  float([-div(sN-1,2):div(sN-1,2)])
yn =  float([-div(sN-1,2):div(sN-1,2)]);
kxmbz =  linspace(-π/q, π/q, sN);

BP.saveplots(ψrex,ψr,ψkex,ψk,ψkmbzex,ψkmbz, xm,yn,float([kx]),float([ky]),kxmbz,q);
BP.matplotspect([freq], ⅈ, spectrum, sω0; vert=true);
BP.matplotcomp(float([ky]),ψkmbz,χ,sζ);
