using PyPlot
import BP

# system parameters
const N = 45
const q = 11
const κ = 0.02
const γ = 0.001
const ν = linspace(-3.45,-2.47,981);
##


# exact spectrum, first 29 eigenvalues
exstates = BP.ExactStates(29, :landau, N, 1/q, κ)


#for plotting filter markers
βlan = [0,1,3,5]
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


δpmp(n₀::Int,m₀::Int) = BP.δpmp(N; n0=n₀, m0=m₀)
gausspmp(n₀::Int,m₀::Int) = BP.gausspmp(N; σ=1., n0=n₀, m0=m₀)
homopmp() = BP.homopmp(N)
randpmp(s::Int) = BP.randpmp(N; seed=s) #1234

           
prm = (1/q,γ,κ);

spδl = BP.Spectrum(ν,δpmp(5,5), :landau, prm...)
spgaussl = BP.Spectrum(ν,gausspmp(5,5), :landau, prm...)
sphoml = BP.Spectrum(ν,homopmp(), :landau, prm...)


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


f, axes = plt.subplots(4, figsize=(10, 7.3))
axes[1][:plot](spδl.νs,spδl.intensity,"k") 
for ω in sω0real
    axes[1][:axvline](x = ω, color="k", ls="dotted")
end
axes[1][:set_ylim](0, 1e3)
axes[1][:yaxis][:set_ticks]([0, 1e3])
axes[1][:yaxis][:set_ticklabels]([L"$0$", L"$10^3$"])

axes[2][:plot](spgaussl.νs,spgaussl.intensity,"k") 
axes[2][:plot](spgausss.νs,spgausss.intensity, color="green", ls="dashed", linewidth=1.5)
for ω in sω0sym
    axes[2][:axvline](x = ω, color="green", ls="dotted")
end 
axes[2][:set_ylim](0, 1e3)
axes[2][:yaxis][:set_ticks]([0, 1e3])
axes[2][:yaxis][:set_ticklabels]([L"$0$", L"$10^3$"])


axes[3][:plot](sphoml.νs,sphoml.intensity,"k") 
axes[3][:plot](sphoms.νs,sphoms.intensity, color="green", ls="dashed", linewidth=1.5)
for ω in sω0lan
    axes[3][:axvline](x = ω, color="k", ls="dotted")
end 
axes[3][:set_ylim](0, 1e6)
axes[3][:yaxis][:set_ticks]([0, 1e6])
axes[3][:yaxis][:set_ticklabels]([L"$0$", L"$10^6$"])


axes[4][:plot](ν,sprandl,"k") 
for ω in exstates.νs
    axes[4][:axvline](x = ω, color="orange", ls="dotted")
end 
axes[4][:set_xlabel](L"$\omega_0 [J]$")
axes[4][:set_ylim](0, 1e6)
axes[4][:yaxis][:set_ticks]([0, 1e6])
axes[4][:yaxis][:set_ticklabels]([L"$0$", L"$10^6$"])

for (i, ax) in enumerate(axes)
    ax[:set_xlim](ν[1], ν[end])

    i != 4 && ax[:set_xticklabels]([])
end 

f[:savefig]("../../figures/selection.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(f)

#TODO: add extra panel for centered gaussian, in linear scale with a zoom

