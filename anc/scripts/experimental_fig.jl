using PyPlot
import BP

#system parameters
sp=1
sq=7
sγ=0.05
sκ=0.2
sN=11
n=5
m=0
sω1=-3.1
sω2=-1.0
dδ=0.001
#spectral range
freq =  sω1:dδ:sω2
#scaling factor for comparing to analytics
sζ=1

#Landau gauge
ft = ( 
    (n,m,a) -> one(Complex{Float64}), (n,m,a) -> one(Complex{Float64}),
    (n,m,a) -> exp(-im*2π*a*m), (n,m,a) -> exp(im*2π*a*m))
ftex =  map(f -> (x, y, z) -> -f(x, y, z), ft);

g =  gcd(sp,sq)
p =  div(sp,g)
q =  div(sq,g)
α =  p/q

#number of nonzero elements in sparse Hamiltonian matrix
nz = BP.countentries(sN)
# we need to convert from sparse to dense in order to use eigvals
M = full(BP.genspmat(ftex...,(n,m,a) -> 1/2*sκ*(n^2+m^2) + zero(Complex{Float64}),sN,nz,α))
#M is the "exact" hamiltonian matrix, without dissipation and pumping
#exact spectrum, first 15 eigenvalues
spectrum =  eigvals(Hermitian(M), 1:15);

#we filter state η
η = 8
β =  η-1
#at energy
sω0=spectrum[η]
#with a δ-like pump
P = BP.δpmp(sN; n0=n, m0=m)

println("[] filtering at energy ω = $(sω0)")

#compute spectrum
function spect()
    v = Array(Float64, length(freq))
    for (j,ω) in enumerate(freq)
        S = BP.genspmat(ft...,(n,m,a) -> ω + im*sγ - 1/2*sκ*(n^2+m^2), sN,nz,α)
        X = S\P
        v[j] = sum(abs2(X))
    end
    v
end
I = spect()

Np =  q*sN+1
kx = linspace(-π, π, Np)
ky = linspace(-π, π, Np)
kxmbz =  linspace(-π/q, π/q, sN)

xm =  float([-div(sN-1,2):div(sN-1,2)])
yn =  float([-div(sN-1,2):div(sN-1,2)])

#no pumping/dissipation
#returns the nth eigenvector
F = reshape(eigfact(Hermitian(M), η:η)[:vectors], sN,sN)
ψrex = abs2(F)
ψkex = abs2(BP.myfft2(F, kx, ky))
ψkmbzex = BP.mbz(ψkex, q,sN)

#matrix of linear system with dissipation
S = BP.genspmat(ft..., (n,m,a) -> sω0 + im*sγ - 1/2*sκ*(n^2+m^2), sN,nz,α)
X = reshape(S\P, sN,sN)
ψr = abs2(X)
ψk = abs2(BP.myfft2(X, kx, ky))
ψkmbz = BP.mbz(ψk, q,sN)

#analytical w.f. in MBZ
function getχ()
    ydata = linspace(-π,π,Np)
    v = Array(Complex{Float64}, Np)
    for (i,y) in enumerate(ydata)
        v[i] = BP.χ(0.,y,α,β)
    end
    radical = sqrt(sqrt(2/q) / (2π*2^β * factorial(β) * 2π*α))
    radical .* v
end

χ = getχ()


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

#plot spectrum
mx = maximum(I)
l=size(ψkmbz)[2]

f, axes = plt.subplots(2,1, figsize=(10, 6))
for (i,ax) in enumerate(axes)
    if i==1
        ax[:plot](freq, I, "k")
        ax[:vlines](spectrum, 0, mx/2, colors="orange", linestyles="dashed")
        ax[:axvline](x = sω0, color="k", ls="dashed")

        ax[:set_xlim](freq[1], freq[end])
        ax[:set_ylim](0, mx)
        ax[:yaxis][:set_ticks]([0.,1.7,3.5])

        ax[:set_xlabel](L"$\omega_0 [J]$")
        ax[:set_ylabel](L"$\sum_{m,n} |a_{m,n}|^2$ [a.u.]")     
    else
        ax[:plot](ky, abs2(χ), "blue", ls="dotted", linewidth=1.5)
        ax[:plot](ky, ψkmbz[:,div(l-1,2)+1]/sζ, "k")
        ax[:plot](ky, ψkmbzex[:,div(l-1,2)+1]/sζ, color="orange", ls="--")
        
        ax[:set_xlim](-π, π)
        ax[:set_xticks]([-π,-π/2,0,π/2,π])
        ax[:set_xticklabels]([L"$-\pi$",L"$-\pi/2$",L"$0$",L"$\pi/2$",L"$\pi$"])        
        ax[:set_ylabel](L"$|\chi_7(0,p_y)|^2$")
        ax[:set_xlabel](L"$p_y$")
        ax[:yaxis][:set_ticks]([0.,0.2,0.4])
    end
end 

f[:savefig]("../../figures/exp_spect.svg", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(f)


#plot w.f. in real and mom space
#TODO: debug assymetry in right column
f, axes = plt.subplots(2,3, figsize=(10, 7.3))
#real space
for (i,ψ) in enumerate((ψrex, ψr))
    ax = axes[i,1]
    im = ax[:imshow](ψ, origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                     extent=[minimum(xm), maximum(xm), minimum(yn), maximum(yn)], aspect=1,
                     vmin=0, vmax=0.1)
    ax[:set_ylabel](L"$n$")
    if i == 2
        ax[:set_xlabel](L"$m$")
        
        cbaxes = f[:add_axes]([0.1, 0.55, 0.22, 0.015])
        cbar = f[:colorbar](im, cax=cbaxes, orientation="horizontal")
        cbar[:set_ticks]([0, 0.05, 0.1])
        cbar[:set_ticklabels]([L"$0$", L"$0.05$", L"$0.1$"])
        cbar[:set_label](L"$|a_{m,n}|^2$", rotation=0, labelpad=-8, y=0.5)
        cbar[:solids][:set_edgecolor]("face")

    else
        ax[:set_xticklabels]([])
    end 
end

#momentum space
for (i,ψ) in enumerate((ψkex, ψk))
    ax = axes[i,2]
    im = ax[:imshow](ψ, origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                     extent=[minimum(kx), maximum(kx), minimum(ky), maximum(ky)], aspect=1,
                     vmin=0, vmax=14)
    ax[:set_ylabel](L"$p_y$")
    if i == 2
        ax[:set_xlabel](L"$p_x$")
        ax[:set_xticks]([-π,0,π])
        ax[:set_xticklabels]([L"$-\pi$",L"$0$",L"$\pi$"])        

        cbaxes = f[:add_axes]([0.42, 0.55, 0.22, 0.015])
        cbar = f[:colorbar](im, cax=cbaxes, orientation="horizontal")
        cbar[:set_ticks]([0, 7, 14])
        cbar[:set_label](L"$|a_{p_x,p_y}|^2$", rotation=0, labelpad=-8, y=0.5)
        cbar[:solids][:set_edgecolor]("face")
    else
        ax[:set_xticklabels]([])
    end
    ax[:set_yticks]([-π,0,π])
    ax[:set_yticklabels]([L"$-\pi$",L"$0$",L"$\pi$"])        
end

#MBZ
for (i,ψ) in enumerate((ψkmbzex, ψkmbz))
    ax = axes[i,3]
    im = ax[:imshow](ψ, origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                     extent=[minimum(kxmbz), maximum(kxmbz), minimum(ky), maximum(ky)], aspect=1/q,
                     vmin=0, vmax=1)
    ax[:set_ylabel](L"$p_y$")
    if i == 2
        ax[:set_xlabel](L"$p_x^0$") 
        ax[:set_xticks]([-π/q,0,π/q])
        ax[:set_xticklabels]([L"$-\pi/q$",L"$0$",L"$\pi/q$"])

        cbaxes = f[:add_axes]([0.74, 0.55, 0.22, 0.015])
        cbar = f[:colorbar](im, cax=cbaxes, orientation="horizontal")
        cbar[:set_ticks]([0, 0.5, 1])
        cbar[:set_ticklabels]([L"$0$", L"$0.5$", L"$1$"])
        cbar[:set_label](L"$|a_{p_x^0,p_y}|^2$", rotation=0, labelpad=-8, y=0.5)
        cbar[:solids][:set_edgecolor]("face")
    else
        ax[:set_xticklabels]([])
    end
    ax[:set_yticks]([-π,0,π])
    ax[:set_yticklabels]([L"$-\pi$",L"$0$",L"$\pi$"])        
end


f[:savefig]("../../figures/experimental.png", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(f)
