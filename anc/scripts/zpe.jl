using HH
using Base.Test


## bw = [bwidth(q) for q in qs]
## bg = [bgap(q) for q in qs]


## hheig = HH.hhladder(α,q)
## ladder = [mean(hheig[:,1]) + (β + 1/2)*κ/(2π*α) for β in 0:14]
## ## using Color
## ## dc = distinguishable_colors(q)
## a = Layer[]
## for i = 1:q
##     push!(a, layer(x=linspace(-π, π, 100), y=hheig[:,i], Geom.line, Theme(default_color=dc[i]))[1] )
## end
## plot(a)
## [ηzpe(q,κ) for q in qs, κ in κs] #


@test_approx_eq_eps(HH.ηzpe(0.02, 11), 0.9368071126414523, 1e-5)
@test_approx_eq_eps(HH.ηzpe(11, 0.02), 0.9368071126414523, 1e-5)


@time HH.ηzpe(0.02, 11);
@time HH.ηzpe(11, 0.02);
