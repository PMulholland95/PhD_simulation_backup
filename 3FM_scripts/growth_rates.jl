using DelimitedFiles, CairoMakie, GLMakie
GLMakie.activate!()

#k(x/y)n = mode number index
kxn = 1
kyn = 1

kx = range(-0.5,0.5, length=21)
ky = range(0.05,0.5, length=10)

#eva = eigenvalues
eva = fm.eigenvalues
#ykx = []

#γkx corresponds to growth(kx, ky=fixed)
#γky corresponds to growth(ky, kx=fixed)

γkx = circshift(map(x->imag(eva[x, kyn][1]), 1:21), (10,0)) 
γky = map(x->imag(eva[kxn, x][1]), 1:10) 

fig = Figure()

fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2]) 

γkxline = lines!(ax1, kx, γkx, color = :red)
γkyline = lines!(ax2, ky, γky, color = :blue)

ax1.title = "γ"
ax1.xlabel = "kx"

ax2.title = "γ"
ax2.xlabel = "ky"

xlims!(ax1, -0.6, 0.6)
xlims!(ax2, 0.0, 0.55)
ylims!(ax1, 0, 0.035)
ylims!(ax2, 0, 0.22)


save("growth_wavenumber.png",fig)

fig

