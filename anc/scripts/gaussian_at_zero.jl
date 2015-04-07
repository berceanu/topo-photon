#S gauge
#gaussian pump @ (0,0)

using PyPlot
import BP

#system parameters
sp=1
sq=11
sγ=0.001
sκ=0.02
sN=45
n=0
m=0
sω1=-3.45
sω2=-2.47
dδ=0.001
freq =  sω1:dδ:sω2

#Symmetric gauge
ft = (
      (n,m,a) -> exp(-im*π*a*n), (n,m,a) -> exp(im*π*a*n),
      (n,m,a) -> exp(-im*π*a*m), (n,m,a) -> exp(im*π*a*m))


g =  gcd(sp,sq)
p =  div(sp,g)
q =  div(sq,g)
α =  p/q

nz = BP.countnonzeros(sN)

pumps = (BP.gausspmp(sN; A=1., σ=1., n0=n, m0=m),
         BP.gausspmp(sN; A=1., σ=5., n0=n, m0=m),
         BP.gausspmp(sN; A=1., σ=10., n0=n, m0=m),
         BP.gausspmp(sN; A=1., σ=20., n0=n, m0=m))

function spect()
    v = Array(Float64, length(pumps),length(freq))
    for (j,ω) in enumerate(freq)
        S = BP.genspmat(ft...,(n,m,a) -> ω + im*sγ - 1/2*sκ*(n^2+m^2), sN,nz,α)
        for (i,P) in enumerate(pumps)
            X = S\P
            v[i,j] = sum(abs2(X))
        end
    end
    v
end 

ⅈ = spect()


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


f, axes = plt.subplots(length(pumps), figsize=(10, 7.3))
for (i, ax) in enumerate(axes)
    I = vec(ⅈ[i,:])
    mx = maximum(I)

    ax[:plot](freq, I, color="green", ls="dotted", linewidth=1.5)
    ax[:set_xlim](freq[1], freq[end])

    if i != length(pumps)
        ax[:set_xticklabels]([])
    else
        ax[:set_xlabel](L"$\omega_0 [J]$")
    end
    ax[:set_yticklabels]([])
#    ax[:set_yticks]([])
#    ax[:set_ylabel](L"$\sum_{m,n} |a_{m,n}|^2$ [a.u.]")
end 

f[:savefig]("../../figures/gaussian_at_zero.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(f)


#full plot range, both in x and y
xm = [-div(sN-1,2):div(sN-1,2)]
yn = [-div(sN-1,2):div(sN-1,2)]


f, axes = plt.subplots(length(pumps), figsize=(5, 10))
for (i, ax) in enumerate(axes)
    data = abs2(pumps[i])
    im = ax[:imshow](reshape(data, sN,sN), origin="upper", ColorMap("gist_heat_r"), interpolation="none",
                    extent=[minimum(xm), maximum(xm), minimum(yn), maximum(yn)],
                     vmin=0, vmax=1)
    ax[:set_xticks]([-20,0,20])
    ax[:set_yticks]([-20,0,20])

    ax[:set_ylabel](L"$n$", labelpad=-10)
    if i == 4 #bottom
        ax[:set_xlabel](L"$m$", labelpad=-6)
    else
        ax[:set_xticklabels]([])
    end

    cbaxes = f[:add_axes]([0.3, 0., 0.4, 0.015])
cbar = f[:colorbar](im, cax=cbaxes, orientation="horizontal")
cbar[:set_ticks]([0, 1])
cbar[:set_ticklabels]([L"$0$", L"$1$"])
cbar[:set_label](L"$|f_{m,n}|^2$", rotation=0, labelpad=-15, y=0.5)
cbar[:solids][:set_edgecolor]("face")

end 


f[:savefig]("../../figures/gaussian_pumps.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(f)

