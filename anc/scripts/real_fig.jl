using PyPlot
import BP

#system parameters
sp=1
sq=11
sγ=0.001
sκ=0.02
sN=45
n=5
m=5


#Landau gauge
ft = ( 
      (n,m,a) -> one(Complex{Float64}), (n,m,a) -> one(Complex{Float64}),
      (n,m,a) -> exp(-im*2π*a*m), (n,m,a) -> exp(im*2π*a*m))
ftex =  map(f -> (x, y, z) -> -f(x, y, z), ft)



g =  gcd(sp,sq)
p =  div(sp,g)
q =  div(sq,g)
α =  p/q

nz = BP.countnonzeros(sN)
# we need to convert from sparse to dense in order to use eigvals
M = full(BP.genspmat(ftex...,(n,m,a) -> 1/2*sκ*(n^2+m^2) + zero(Complex{Float64}),sN,nz,α))
#M is the "exact" hamiltonian matrix, without dissipation and pumping
#exact spectrum
spectrum =  eigvals(Hermitian(M), 1:29)


P = BP.δpmp(sN; A=1., n0=n, m0=m)


#filtering @
β = [0,6,15,26]
#we filter state η
η = β + 1
#at energy
sω0= [spectrum[state]::Float64 for state in η]



function wfreal(ω)
    #matrix of linear system with dissipation
    S = BP.genspmat(ft..., (n,m,a) -> ω + im*sγ - 1/2*sκ*(n^2+m^2), sN,nz,α)
    X = reshape(S\P, sN,sN)
    abs2(X)
end 

#calculating all real space wfs
ψ = Array(Float64, (sN, sN, length(β)))
for (i,ω) in enumerate(sω0)
    ψ[:,:,i] = wfreal(ω)
end





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

#full plot range, both in x and y
xm = [-div(sN-1,2):div(sN-1,2)]

#zoom in :)
edge = 10
st = findin(xm, -edge)[1]
en = findin(xm, edge)[1]


#plot w.f. in real space
f, axes = plt.subplots(1,length(β), figsize=(10, 5))

for (i,ax) in enumerate(axes)
    im = ax[:imshow](ψ[st:en,st:en,i], origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                     extent=[-edge, edge, -edge, edge])
    ax[:set_xlabel](L"$m$")
    if i == 1 #leftmost panel
        ax[:set_ylabel](L"$n$")
    else
        ax[:set_yticklabels]([])
    end
end 


f[:savefig]("../../figures/real.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(f)

for i = 1:4
    println(maximum(ψ[:,:,i]))
end 
