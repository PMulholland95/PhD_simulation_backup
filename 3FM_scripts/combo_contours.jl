using DelimitedFiles, CairoMakie, GLMakie
GLMakie.activate!()

kx = range(-0.5,0.5, length = 21)
ky = range(0.05,0.5, length = 10)

kxn = LinRange(1,21,21)
kyn = LinRange(1,10,10)

eva = fm.eigenvalues


omg = circshift([real(eva[x,y][1]) for x in 1:21, y in 1:10], (10,0)) 
gam = circshift([imag(eva[x,y][1]) for x in 1:21, y in 1:10], (10,0)) 
marg = circshift([real(eva[x,y][3]) for x in 1:21, y in 1:10], (10,0)) 


f = Figure()

fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

Axis(f[1, 1], title = "ωᵣ - (un)stable modes", xlabel = "kx", ylabel = "ky", xticks = -0.5:0.5:11)
cofreq = contourf!(kx, ky, omg, levels = 25, colormap = :vik)

Axis(f[1, 3], title = "γ - unstable modes", xlabel = "kx", ylabel = "ky", xticks = -0.5:0.5:11)
cogam = contourf!(kx, ky, gam, levels = 25, colormap = :vik)

Axis(f[2, 1], title = "ωᵣ - marginal modes", xlabel = "kx", ylabel = "ky", xticks = -0.5:0.5:11)
comarg = contourf!(kx, ky, marg, levels = 25, colormap = :vik)

Colorbar(f[1, 2], cofreq)

Colorbar(f[1, 4], cogam)

Colorbar(f[2, 2], comarg)


save("all_contour.png",f)

f


