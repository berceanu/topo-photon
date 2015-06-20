using PyPlot
using PyCall
@pyimport mpl_toolkits.axes_grid1.inset_locator as axloc

import BP

# system parameters
const N = 45
const q = 11
const κ = 0.02
const γ = 0.001
const ν = linspace(-3.45,-2.47,981)
##


# exact spectrum, first 29 eigenvalues
exstates = BP.ExactStates(29, :landau, N, 1/q, κ)


#for plotting filter markers
βlan = [0,2,4,6]
βsym = [0,1,9,20]
βreal = [0,6,15,26]

#we filter state η
ηlan = βlan + 1
ηsym = βsym + 1
ηreal = βreal + 1

#at energy
sω0lan = [exstates.νs[state]::Float64 for state in ηlan]
sω0sym = [exstates.νs[state]::Float64 for state in ηsym]
sω0real= [exstates.νs[state]::Float64 for state in ηreal]

# energy boundaries for state with β=4
ω₁= exstates.νs[4] - 0.005
ω₂= exstates.νs[6] + 0.005


δpmp(n₀::Int,m₀::Int) = BP.δpmp(N; n0=n₀, m0=m₀)
gausspmp(n₀::Int,m₀::Int) = BP.gausspmp(N; σ=1., n0=n₀, m0=m₀)
homopmp() = BP.homopmp(N)
randpmp(s::Int) = BP.randpmp(N; seed=s) #1234

prm = (1/q,γ,κ);

spδl = BP.Spectrum(ν,δpmp(5,5), :landau, prm...)
spgaussl = BP.Spectrum(ν,gausspmp(5,5), :landau, prm...)
sphoml = BP.Spectrum(ν,homopmp(), :landau, prm...)

## deprecated ##
# centered gaussian
# spgauss0 = BP.Spectrum(ν, BP.gausspmp(N; A=1., σ=20., n0=0, m0=0), :symmetric, prm...)


spgausss = BP.Spectrum(ν,gausspmp(5,5), :symmetric, prm...)
sphoms = BP.Spectrum(ν,homopmp(), :symmetric, prm...)

## averaging over 100 random phase distributions ##
# reading result from file | it takes a while to compute ;)
sprandl = vec(readdlm("sprandl.txt", Float64))

# computing result
## intvec = zeros(Float64, length(ν));
## A = spzeros(Complex{Float64}, N^2,N^2);
## for j=1:100
##     P=randpmp(j)
##     for (i,ω) in enumerate(ν)
##         BP.buildham_landau!(A, N,1/q,κ,γ, ω)
##         intvec[i] += sum(abs2(A\P))
##     end 
## end 
## sprandl = intvec./100;

#extrema(sprandl)

# writing computed result to file
## writedlm("sprandl.txt", sprandl)


##
for sp in (spδl,spgaussl,sphoml)
    println(extrema(sp.intensity))
end

println()

for sp in (spgausss, sphoms)
    println(extrema(sp.intensity))
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


# half the level spacing
hf = (exstates.νs[2] - exstates.νs[1])/2

fig, axes = plt.subplots(4, figsize=(10, 7.3))
axes[1][:plot](spδl.νs,spδl.intensity,"k") 
for (i,ω) in enumerate(sω0real)
    axes[1][:axvline](x = ω, color="orange", ls="dotted")
    axes[1][:text](ω + hf/4, 8e2, string(βreal[i]))
end
axes[1][:set_ylim](0, 1e3)
axes[1][:yaxis][:set_ticks]([0, 1e3])
axes[1][:yaxis][:set_ticklabels]([L"$0$", L"$10^3$"])
axes[1][:text](ν[1] + hf, 5e2, "(a)")


axes[2][:plot](spgaussl.νs,spgaussl.intensity,"k") 
axes[2][:plot](spgausss.νs,spgausss.intensity, color="green", ls="dashed", linewidth=1.5)
for (i,ω) in enumerate(sω0sym)
    axes[2][:axvline](x = ω, color="orange", ls="dotted")
    axes[2][:text](ω + hf/4, 8e2, string(βsym[i]))
end 
axes[2][:set_ylim](0, 1e3)
axes[2][:yaxis][:set_ticks]([0, 1e3])
axes[2][:yaxis][:set_ticklabels]([L"$0$", L"$10^3$"])
axes[2][:text](ν[1] + hf, 5e2, "(b)")


axes[3][:plot](sphoml.νs,sphoml.intensity,"k") 
axes[3][:plot](sphoms.νs,sphoms.intensity, color="green", ls="dashed", linewidth=1.5)
# insert with zoom of peak β=4
axins = axloc.inset_axes(axes[3],
                        width="30%", # width = 30% of parent_bbox
                        height="50%",
                        loc=9) # located at upper middle part
# plot same thing as in parent box
axins[:plot](sphoms.νs,sphoms.intensity, color="green", ls="dashed", linewidth=1.5)
axins[:plot](sphoml.νs,sphoml.intensity, color="black") 

# but set much narrower limits
axins[:set_xlim]([ω₁, ω₂])
axins[:xaxis][:set_ticks]([ω₁, ω₂])
axins[:xaxis][:set_ticklabels]([string(round(ω₁,2)), string(round(ω₂,2))], fontsize=8)


axins[:set_ylim]([1000, 10000])
axins[:yaxis][:set_ticks]([1000, 10000])
axins[:yaxis][:set_ticklabels]([L"$10^3$", L"$10^4$"], fontsize=8)

# draw vertical lines at position of every exact eigenstate
for (i,ω) in enumerate(exstates.νs[4:6])
    axins[:axvline](x = ω, color="orange", ls="dotted")
    axins[:text](ω + hf/8, 8e3, string(i+2), fontsize=8)
end 

# transparency setting, not needed when exporting to pdf
#axins[:patch][:set_alpha](1.0)

# hide inset tick labels
#plt.xticks(visible=false)
#plt.yticks(visible=false)

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
axloc.mark_inset(axes[3], axins, loc1=2, loc2=4, ec="0.", fc="none")

for (i,ω) in enumerate(sω0lan)
    axes[3][:axvline](x = ω, color="orange", ls="dotted")
    axes[3][:text](ω + hf/4, 8e5, string(βlan[i]))
end 
axes[3][:set_ylim](0, 1e6)
axes[3][:yaxis][:set_ticks]([0, 1e6])
axes[3][:yaxis][:set_ticklabels]([L"$0$", L"$10^6$"])
axes[3][:text](ν[1] + hf, 5e5, "(c)")


axes[4][:plot](ν,sprandl,"k") 
for ω in exstates.νs
    axes[4][:axvline](x = ω, color="orange", ls="dotted")
end 
axes[4][:set_xlabel](L"$\omega_0 [J]$")
axes[4][:set_ylim](0, 1e6)
axes[4][:yaxis][:set_ticks]([0, 1e6])
axes[4][:yaxis][:set_ticklabels]([L"$0$", L"$10^6$"])
axes[4][:text](ν[1] + hf, 5e5, "(d)")

for (i, ax) in enumerate(axes)
    ax[:set_xlim](ν[1], ν[end])
    i != 4 && ax[:set_xticklabels]([])
end 

# set common y label to all subplots
fig[:text](0.0, 0.5, L"$\sum_{m,n} |a_{m,n}|^2$ [a.u.]", ha="center", va="center", rotation="vertical")

fig[:savefig]("../../figures/selection.pdf", pad_inches=0.0, bbox_inches="tight")
plt.close(fig)



## deprecated ##
# plotting panel for centered gaussian
## fig, ax = plt.subplots(figsize=(8, 3))
##     ax[:plot](spgauss0.νs, spgauss0.intensity, color="green", linewidth=1.5)

##     ax[:set_xlim](spgauss0.νs[1], spgauss0.νs[end])

##     ax[:set_xlabel](L"$\omega_0 [J]$")
##     ax[:set_yticklabels]([])

## fig[:savefig]("../../figures/gaussian.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
## plt.close(fig)

# TODO: merge with contents of real_fig.jl
# TODO: merge with contents of momentum_fig.jl
