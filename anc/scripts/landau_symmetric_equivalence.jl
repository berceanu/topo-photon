#both gauges
#fixed β
#δ pump @ (5,5)
#L: R I ||^2
#S: R I ||^2

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
ftlan = ( 
         (n,m,a) -> one(Complex{Float64}), (n,m,a) -> one(Complex{Float64}),
         (n,m,a) -> exp(-im*2π*a*m), (n,m,a) -> exp(im*2π*a*m))
ftex =  map(f -> (x, y, z) -> -f(x, y, z), ftlan)


#Symmetric gauge
ftsym = (
         (n,m,a) -> exp(-im*π*a*n), (n,m,a) -> exp(im*π*a*n),
         (n,m,a) -> exp(-im*π*a*m), (n,m,a) -> exp(im*π*a*m))

g =  gcd(sp,sq)
p =  div(sp,g)
q =  div(sq,g)
α =  p/q

nz = BP.countentries(sN)
# we need to convert from sparse to dense in order to use eigvals
M = full(BP.genspmat(ftex...,(n,m,a) -> 1/2*sκ*(n^2+m^2) + zero(Complex{Float64}),sN,nz,α))
#M is the "exact" hamiltonian matrix, without dissipation and pumping
#exact spectrum
spectrum =  eigvals(Hermitian(M), 1:29)


P = BP.δpmp(sN; A=1., n0=n, m0=m)

β = 15
#we filter state η
η = β + 1
#at energy
sω0 = spectrum[η]



function wfreal(gauge)
    #matrix of linear system with dissipation
    S = BP.genspmat(gauge..., (n,m,a) -> sω0 + im*sγ - 1/2*sκ*(n^2+m^2), sN,nz,α)
    reshape(S\P, sN,sN)
end 


ψ = Array(Float64, (sN, sN, 2, 2))

for (i,ft) in enumerate((ftlan,ftsym))
    X = wfreal(ft)
    ψ[:,:,1,i] = abs(X)
    ψ[:,:,2,i] = angle(X)
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



f, axes = plt[:subplots](2,2, figsize=(5,4.7))
for i = 1:2 #loop over columns
    #top row: L
    ax = axes[1,i]
    im = ax[:imshow](ψ[st:en,st:en,i,1], origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                     extent=[-edge,edge,-edge,edge])
    ax[:set_xticklabels]([])
    ax[:set_xticks]([-edge,0,edge])
    ax[:set_yticks]([-edge,0,edge])
    if i == 1 #leftmost panel
        ax[:set_ylabel](L"$n$")
    else
        ax[:set_yticklabels]([])
    end

    #bottom row: S
    ax = axes[2,i]
    im = ax[:imshow](ψ[st:en,st:en,i,2], origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                     extent=[-edge,edge,-edge,edge])
    ax[:set_xlabel](L"$m$")
    ax[:set_xticks]([-edge,0,edge])
    ax[:set_yticks]([-edge,0,edge])
    if i == 1 #leftmost panel
        ax[:set_ylabel](L"$n$")
    else
        ax[:set_yticklabels]([])
    end
end 

f[:tight_layout]()
f[:savefig]("../../figures/equivalence.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt[:close](f)

#in absolute value they are the same, but the phase is different
using Base.Test
@test_approx_eq_eps(ψ[:,:,1,1], ψ[:,:,1,2], 1e-10)

#ranges: abs|0,10 phase|-π,π




