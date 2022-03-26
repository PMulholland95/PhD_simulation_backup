using DelimitedFiles, CairoMakie, GLMakie
GLMakie.activate!()

sc = readdlm("p04sc1_scan_kjm20_beta_1.1%")

markersize = 50
linewidth = 20
fontsize = 80
resolution = (2000,1000)

gam = sc[:,5];
gam = convert(Array{Float64}, gam);

omg = sc[:,6];
omg = convert(Array{Float64}, omg);

ky = sc[:,3];
ky = convert(Array{Float64}, ky);

xax = 0:0.05:0.25
xax = convert(Array{Float64}, xax);
yax = 0:0.02:0.12
yax = convert(Array{Float64}, yax);
yax2 = 0:0.005:0.03
yax2 = convert(Array{Float64}, yax2);
	    
fig = Figure()

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

gl = scatter!(ax2, ky, gam, color = :red, markersize=markersize, marker = '□')
ol = scatter!(ax1, ky, omg, color = :blue, markersize=markersize
, marker = '□')

hlines!(ax1, [0], color = :black, linewidth = linewidth) 
vlines!(ax1, [0], color = :black, linewidth = linewidth) 
hlines!(ax1, yax, color = :grey)
vlines!(ax1, xax, color = :grey)

hlines!(ax2, [0], color = :black, linewidth = linewidth) 
vlines!(ax2, [0], color = :black, linewidth = linewidth)
hlines!(ax2, yax2, color = :grey)
vlines!(ax2, xax, color = :grey)

ax1.title = "β=1.1%, omtᵢ=1.6, omtₑ=omn=0"
ax2.ylabel = "γ" 
ax1.ylabel = "ωᵣ"

ax2.xlabel = "kyρ"
ax1.xlabel = "kyρ"

fig
