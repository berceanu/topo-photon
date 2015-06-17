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

#BP.δpmp, BP.gausspmp, BP.homopmp, BP.randpmp
#(sN; A=1., seed=1234, σ=tbσ, n0=n, m0=m)
Plan = BP.homopmp(sN)
#Psym = BP.homopmp(sN)
Psym = BP.gausspmp(sN; A=1., σ=tbσ, n0=n, m0=m)


βlan = [0,2,4,6]
βsym = [0,1,9,20]
#we filter state η
ηlan = βlan + 1
ηsym = βsym + 1
#at energy
sω0lan = [spectrum[state]::Float64 for state in ηlan]
sω0sym = [spectrum[state]::Float64 for state in ηsym]


Np = 200 #no MBZ here
kx = linspace(-π, π, Np)
ky = linspace(-π, π, Np)


function wfmom(gauge, ω, P)
    #matrix of linear system with dissipation
    S = BP.genspmat(gauge..., (n,m,a) -> ω + im*sγ - 1/2*sκ*(n^2+m^2), sN,nz,α)
    X = reshape(S\P, sN,sN)
    abs2(BP.myfft2(X, kx, ky))
end 

#calculating all mom space wfs
ψL = Array(Float64, (Np, Np, length(βlan)))
ψS = Array(Float64, (Np, Np, length(βsym)))
for (i,ω) in enumerate(sω0lan)
    ψL[:,:,i] = wfmom(ftlan, ω, Plan)
end
for (i,ω) in enumerate(sω0sym)
    ψS[:,:,i] = wfmom(ftsym, ω, Psym)
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

#plot w.f. in  mom space

f, axes = plt.subplots(2,length(βlan), figsize=(10, 5))

for i = 1:length(βlan) #loop over columns
    #top row
    ax = axes[1,i]
    im = ax[:imshow](ψL[:,:,i], origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                     extent=[minimum(kx), maximum(kx), minimum(ky), maximum(ky)])
    ax[:set_xticklabels]([])
    ax[:set_xticks]([-π,0,π])
    ax[:set_yticks]([-π,0,π])
    if i == 1 #leftmost panel
        ax[:set_ylabel](L"$p_y$")
        ax[:set_yticklabels]([L"$-\pi$",L"$0$",L"$\pi$"])
    else
        ax[:set_yticklabels]([])
    end

    #bottom row
    ax = axes[2,i]
    if i == 3 #third pannel
        im = ax[:imshow](ψS[:,:,i], origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                         extent=[minimum(kx), maximum(kx), minimum(ky), maximum(ky)],
                         vmin=0, vmax=270000)
        
    else
        im = ax[:imshow](ψS[:,:,i], origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                         extent=[minimum(kx), maximum(kx), minimum(ky), maximum(ky)])
    end
    ax[:set_xlabel](L"$p_x$")
    ax[:set_xticks]([-π,0,π])
    ax[:set_yticks]([-π,0,π])
    ax[:set_xticklabels]([L"$-\pi$",L"$0$",L"$\pi$"])
    if i == 1 #leftmost panel
        ax[:set_ylabel](L"$p_y$")
        ax[:set_yticklabels]([L"$-\pi$",L"$0$",L"$\pi$"])
    else
        ax[:set_yticklabels]([])
    end
end 


f[:savefig]("../../figures/momentum.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(f)



println("Landau g. extrema (top row):")
for i = 1:length(βlan) #loop over columns
    println(extrema(ψL[:,:,i]))
end

println("Symmetric g. extrema (bottom row):")
for i = 1:length(βlan) #loop over columns
    println(extrema(ψS[:,:,i]))
end


