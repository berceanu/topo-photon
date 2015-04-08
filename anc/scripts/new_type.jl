import BP
#reload("BP")

import PyPlot; const plt = PyPlot
using LaTeXStrings



δpmp(n₀::Int,m₀::Int) = BP.δpmp(45; n0=n₀, m0=m₀)

P = δpmp(5,5)

@time a = BP.WaveFunction(-2.0,P,:symmetric)

b = BP.WaveFunction(-2.0,P)

## sω1=-2.49
## sω2=-2.47
## dδ=0.001
## freq =  [sω1:dδ:sω2]

freq = linspace(-3,-2,1000)


@time sp = BP.Spectrum(freq,P);

BP.getstate(sp, -2.485)

# matplotlib parameters
plt.matplotlib["rcParams"][:update](["axes.labelsize" => 22,
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
f, ax = plt.subplots(figsize=(10, 4))

ax[:plot](sp.νs, sp.intensity, "k")
ax[:axvline](x = -2.5, color="k", ls="dashed")

ax[:set_xlim](sp.νs[1], sp.νs[end])

ax[:set_xlabel](L"$\omega_0 [J]$")
ax[:set_ylabel](L"$\sum_{m,n} |a_{m,n}|^2$ [a.u.]")     

f[:savefig]("spect.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(f)



## X = reshape(S\P, N,N)
## k = linspace(-π, π, 200)
## Xk = myfft2(X, k,k)
## abs(), abs2(), angle()










    
