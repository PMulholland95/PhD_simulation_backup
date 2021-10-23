using DelimitedFiles, CairoMakie, GLMakie
GLMakie.activate!()

kx = range(-0.5,0.5, length = 21)
ky = range(0.05,0.5, length = 10)

kxn = LinRange(1,21,21)
kyn = LinRange(1,10,10)

eva = fm.eigenvalues

gam = circshift([imag(eva[x,y][1]) for x in 1:21, y in 1:10], (10,0)) 

f = Figure()

fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

Axis(f[1, 1], title = "γ", xlabel = "kx", ylabel = "ky", xticks = -0.5:0.5:11)
co = contourf!(kx, ky, gam, levels = 25, colormap = :vik)

Colorbar(f[1, 2], co)

save("growth_wavenumber_contour.png",f)

f


