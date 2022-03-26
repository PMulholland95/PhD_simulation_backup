using DelimitedFiles, CairoMakie, LaTeXStrings 

markersize = 80
linewidth = 15
fontsize = 120
resolution = (2000,2000)

fig = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax = Axis(fig[1,1])

xs = range(0, 10, length = 30)

points = [Point2f(x, y) for y in -10:1:10 for x in -10:1:10]
rotations = range(0, 2pi, length = length(points))

scatter!(points, rotations = rotations, markersize = markersize, marker = '↑')

save("flow_test.png",fig)

fig
