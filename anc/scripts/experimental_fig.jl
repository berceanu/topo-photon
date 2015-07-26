import BP
using Base.Test
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

# system parameters
sN=11 # true system size

# where we pump
n=5
m=0

# inverse lifetime
sγ=0.05

# strength of the trap
sκ = 0.2

q = 7
r = 11 # points in MBZ

# N should be an odd multiple of q
N = r*q # zero-padded system size
l = div(N-1,2)
x = -l:l
δk = 2π/N #resolution in mom space
k = x * δk

@test x[l+1] == 0
@test k[l+1] == 0.0
@test isodd(N)
@test isodd(r)

#generate all p values inside MBZ
xmbz = -div(r-1,2):div(r-1,2)
kxmbz = xmbz * δk

# no pumping or dissipation

# exact spectrum, first 15 eigenvalues
exexp = BP.ExactStates(15, :landau, sN, 1/q, sκ)


#we filter state η
η = 8
β =  η-1
#at energy
sω0=exexp.νs[η]
println("[ ] filtering at energy ω = $(sω0)")

#with a δ-like pump
δpmp(n₀::Int,m₀::Int) = BP.δpmp(sN; n0=n₀, m0=m₀)
P = δpmp(n,m)

# compute spectrum
sω1=-3.1
sω2=-1.0
dδ=0.001
#spectral range
ν =  sω1:dδ:sω2

sp = BP.Spectrum([ν], P, :landau, 1/q, sγ, sκ)


# no pumping or dissipation
state = BP.getstate(exexp, η) #11x11 Array{Complex{Float64},2} (0.0, 0.025)
ψrex = abs2(state)
ψkex = abs2(BP.myfft2(state, k, k)) #77x77 Array{Float64,2} (0.0, 21.3)
ψkmbzex = BP.mbz(ψkex, r, q, kxmbz, k) #77x11 Array{Float64,2} (0., 1.24)

# with pumping
X = BP.getstate(sp, sω0)
ψr = abs2(X)
ψk = abs2(BP.myfft2(X, k, k))
ψkmbz = BP.mbz(ψk, r, q, kxmbz, k)

#analytical w.f. in MBZ
function getχ(Np,q,β)
    α = 1/q
    ydata = linspace(-π,π,Np)
    v = Array(Complex{Float64}, Np)
    for (i,y) in enumerate(ydata)
        v[i] = BP.χ(0.,y,α,β)
    end
    radical = sqrt(sqrt(2/q) / (2π*2^β * factorial(β) * 2π*α))
    radical .* v
end

χ = getχ(N,q,β)

# plotting

#plot spectrum
mx = maximum(sp.intensity)

fig, axes = plt.subplots(2,1, figsize=(10, 6))
for (i,ax) in enumerate(axes)
    if i==1
        ax[:plot](sp.νs, sp.intensity, "k")
        ax[:vlines](exexp.νs, 0, mx/2, colors="orange", linestyles="dashed")
        ax[:axvline](x = sω0, color="k", ls="dashed")

        ax[:set_xlim](sp.νs[1], sp.νs[end])
        ax[:set_ylim](0, mx)
        ax[:yaxis][:set_ticks]([0.,1.7,3.5])

        ax[:set_xlabel](L"$\omega_0/J$")
        ax[:set_ylabel](L"$\sum_{m,n} |a_{m,n}|^2$ [a.u.]")     
    else
        ax[:plot](k, abs2(χ), "blue", ls="dotted", linewidth=1.5)
        ax[:plot](k, ψkmbz[:,6], "k")
        ax[:plot](k, ψkmbzex[:,6], color="orange", ls="--")
        
        ax[:set_xlim](-π, π)
        ax[:set_xticks]([-π,-π/2,0,π/2,π])
        ax[:set_xticklabels]([L"$-\pi$",L"$-\pi/2$",L"$0$",L"$\pi/2$",L"$\pi$"])        
        ax[:set_ylabel](L"$|\chi_7(0,p_y)|^2$")
        ax[:set_xlabel](L"$p_y$")
        ax[:yaxis][:set_ticks]([0.,0.2,0.4])
    end
end 

fig[:savefig]("../../figures/exp_spect.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(fig)


#plot w.f. in real and mom space
fig, axes = plt.subplots(2,3, figsize=(10, 7.3))
#real space
for (i,ψ) in enumerate((ψrex, ψr))
    ax = axes[i,1]
    img = ax[:imshow](ψ, origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                     extent=[-5.5, 5.5, -5.5, 5.5], aspect=1,
                     vmin=0, vmax=0.1)
    ax[:set_ylabel](L"$n$")
    if i == 2
        ax[:set_xlabel](L"$m$")
        
        cbaxes = fig[:add_axes]([0.1, 0.55, 0.22, 0.015])
        cbar = fig[:colorbar](img, cax=cbaxes, orientation="horizontal")
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
    img = ax[:imshow](ψ, origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                     extent=[-π, π, -π, π], aspect=1,
                     vmin=0, vmax=14)
    ax[:set_ylabel](L"$p_y$")
    ax[:set_xticks]([-π,0,π])
    
    if i == 2
        ax[:set_xlabel](L"$p_x$")
        ax[:set_xticklabels]([L"$-\pi$",L"$0$",L"$\pi$"])        

        cbaxes = fig[:add_axes]([0.42, 0.55, 0.22, 0.015])
        cbar = fig[:colorbar](img, cax=cbaxes, orientation="horizontal")
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
    img = ax[:imshow](ψ, origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                     extent=[-π/q, π/q, -π, π], aspect=1/q,
                     vmin=0, vmax=1)
    ax[:set_ylabel](L"$p_y$")
    ax[:set_xticks]([-π/q,0,π/q])
    
    if i == 2
        ax[:set_xlabel](L"$p_x^0$") 
        ax[:set_xticklabels]([L"$-\pi/7$",L"$0$",L"$\pi/7$"])

        cbaxes = fig[:add_axes]([0.74, 0.55, 0.22, 0.015])
        cbar = fig[:colorbar](img, cax=cbaxes, orientation="horizontal")
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


fig[:savefig]("../../figures/experimental.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(fig)

# for debugging

## fig, ax = plt.subplots(figsize=(5, 5))
## ax[:imshow](ψrex, origin="upper", ColorMap("gist_heat_r"), interpolation="none",
##                  extent=[-5.5, 5.5, -5.5, 5.5],
##                  vmin=0, vmax=0.1)

## ax[:set_ylabel](L"$n$", labelpad=-10)
## ax[:set_xlabel](L"$m$", labelpad=-6)
    

## fig, ax = plt.subplots(figsize=(5, 5))
## ax[:imshow](ψkex,
##             origin="upper", ColorMap("gist_heat_r"), interpolation="none",
##             extent=[minimum(k), maximum(k), minimum(k), maximum(k)],
##             vmin=0, vmax=14)

## ax[:set_xticks]([-π,0,π])
## ax[:set_yticks]([-π,0,π])
## ax[:set_xticklabels]([L"$-\pi$",L"$0$",L"$\pi$"])
## ax[:set_yticklabels]([L"$-\pi$",L"$0$",L"$\pi$"])
## ax[:set_xlabel](L"$p_x$", labelpad=-4)
## ax[:set_ylabel](L"$p_y$", labelpad=-10)

##


## fig, ax = plt.subplots(figsize=(5, 5))

## ax[:imshow](ψkmbzex, origin="upper", ColorMap("gist_heat_r"), interpolation="none",
##             extent=[minimum(kxmbz), maximum(kxmbz), minimum(k), maximum(k)], aspect=1/7,
##             vmin=0, vmax=1)
            
## ax[:set_xlabel](L"$p_x^0$") 
## ax[:set_xticks]([-π/7,0,π/7])
## ax[:set_xticklabels]([L"$-\pi/7$",L"$0$",L"$\pi/7$"])

## ax[:set_ylabel](L"$p_y$")
## ax[:set_yticks]([-π,0,π])
## ax[:set_yticklabels]([L"$-\pi$",L"$0$",L"$\pi$"])


#TODO: plot seccond ladder vertical lines with different color
#TODO: plot vertical lines for y slices
#TODO: combine panels in unique figure
