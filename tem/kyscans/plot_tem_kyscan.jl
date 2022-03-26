using DelimitedFiles, CairoMakie 

sc = readdlm("combo_kyscan")

markersize = 80
linewidth = 15
fontsize = 120
resolution = (4000,2200)

gam = sc[:,5];
gam = convert(Array{Float64}, gam);

omg = sc[:,6];
omg = convert(Array{Float64}, omg);

ky = sc[:,3];
ky = convert(Array{Float64}, ky);
	    
fig = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], xlabel="kyρ", ylabel="ωᵣ", linewidth=linewidth)
ax2 = Axis(fig[1,2], xlabel="kyρ", ylabel="γ", linewidth=linewidth)  

lines!(ax2, ky, gam, color = :red,linewidth=linewidth)
lines!(ax1, ky, omg, color = :blue, linewidth=linewidth)

hlines!(ax1, [0], color = :black, linewidth = linewidth) 
vlines!(ax1, [0], color = :black, linewidth = linewidth)
hlines!(ax2, [0], color = :black, linewidth = linewidth) 
vlines!(ax2, [0], color = :black, linewidth = linewidth)

save("tem_kyscan.png", fig)

fig
