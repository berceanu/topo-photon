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
tbσ=1.0
sω1=-3.45
sω2=-2.47
dδ=0.001
#Landau gauge
ftlan = ( 
         (n,m,a) -> one(Complex{Float64}), (n,m,a) -> one(Complex{Float64}),
         (n,m,a) -> exp(-im*2π*a*m), (n,m,a) -> exp(im*2π*a*m))


freq =  sω1:dδ:sω2
ftex =  map(f -> (x, y, z) -> -f(x, y, z), ftlan)

# Is `spectrum` affected by the gauge? NO
# Is symmetric gauge spectrum with a gaussian *identical* to landau gauge with delta? NO

#Symmetric gauge
ftsym = (
         (n,m,a) -> exp(-im*π*a*n), (n,m,a) -> exp(im*π*a*n),
         (n,m,a) -> exp(-im*π*a*m), (n,m,a) -> exp(im*π*a*m))

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

pumps = (BP.δpmp, BP.gausspmp, BP.homopmp, BP.randpmp)

function spect()
    vlan = Array(Float64, length(pumps),length(freq))
    vsym = Array(Float64, length(pumps),length(freq))
    for (j,ω) in enumerate(freq)
        Slan = BP.genspmat(ftlan...,(n,m,a) -> ω + im*sγ - 1/2*sκ*(n^2+m^2), sN,nz,α)
        Ssym = BP.genspmat(ftsym...,(n,m,a) -> ω + im*sγ - 1/2*sκ*(n^2+m^2), sN,nz,α)
        for (i,pump) in enumerate(pumps)
            P = pump(sN; A=1., seed=1234, σ=tbσ, n0=n, m0=m)
            Xlan = Slan\P
            Xsym = Ssym\P
            vlan[i,j] = sum(abs2(Xlan))
            vsym[i,j] = sum(abs2(Xsym))
        end
    end
    (vlan, vsym)
end 

(ⅈlan, ⅈsym) = spect()


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


f, axes = plt.subplots(length(pumps),1, figsize=(10, 7.3))

for (i, ax) in enumerate(axes)
    Ilan = vec(ⅈlan[i,:])
    Isym = vec(ⅈsym[i,:])
    mx = maximum(Ilan)

    ax[:plot](freq, Ilan, "k")
    ax[:plot](freq, Isym, color="green", ls="dotted", linewidth=1.5)

    
    ax[:set_xlim](freq[1], freq[end])
    ax[:set_ylim](0, mx/2)

    if i != length(pumps)
        ax[:set_xticklabels]([])
    else
        ax[:vlines](spectrum, 0, mx/2, colors="orange", linestyles="dashed")
        ax[:set_xlabel](L"$\omega_0 [J]$")
    end
    ax[:set_yticklabels]([])
    ax[:set_yticks]([])
#   ax[:set_ylabel](L"$\sum_{m,n} |a_{m,n}|^2$ [a.u.]")
end 

f[:savefig]("../../figures/selection.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(f)

#TODO: use subplots_adjust() to add space to the left side of the figure
#TODO: use figtext() to add common y label
#TODO: add labels for 0 and maximum intensity
#TODO: add (a), (b), (c), (d)
