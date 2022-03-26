using CairoMakie

resolution=(4000,600)
colormap=:vik
clevels=20

xs = LinRange(0, 10, 500)
ys = LinRange(0, 10, 500)
zs = [cos(2x) * y for x in xs, y in ys]

f = Figure(resolution=resolution)
Axis(f[1, 1])

co = contourf!(xs, ys, zs, levels = clevels, colormap = colormap)

Colorbar(f[1, 2], co)

f
	
