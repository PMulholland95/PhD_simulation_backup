using DelimitedFiles, CairoMakie, GLMakie
GLMakie.activate!()

kx = range(-0.5,0.5, length = 21)
ky = range(0.05,0.5, length = 10)

eva = fm.eigenvalues

marg = circshift([real(eva[x,y][3]) for x in 1:21, y in 1:10], (10,0)) 

f = Figure()

fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

Axis(f[1, 1], title = "ωᵣ - marginal modes", xlabel = "kx", ylabel = "ky")
co = contourf!(kx, ky, marg, levels = 25, colormap = :vik)

Colorbar(f[1, 2], co)

save("marg_freq_contour.png",f)

f


