using DelimitedFiles, CairoMakie 

sc = readdlm("omn_ky_scan")

markersize = 80
linewidth = 5
fontsize = 20
resolution = (4000,2200)


n=23

gam = sc[:,7];
gam = convert(Array{Float64}, gam);

omg = sc[:,8];
omg = convert(Array{Float64}, omg);

omn = sc[1:23,3];
omn = convert(Array{Float64}, omn); 

g1 = gam[1:n];
g2 = gam[1+n:2n];
g3 = gam[1+2n:3n];
g4 = gam[1+3n:4n];
g5 = gam[1+4n:5n];
g6 = gam[1+5n:6n];
g7 = gam[1+6n:7n];

o1 = omg[1:n];
o2 = omg[1+n:2n];
o3 = omg[1+2n:3n];
o4 = omg[1+3n:4n];
o5 = omg[1+4n:5n];
o6 = omg[1+5n:6n];
o7 = omg[1+6n:7n];

gmat = [g1 g2 g3 g4 g5 g6 g7];
omat = [o1 o2 o3 o4 o5 o6 o7];

ky = 0.4,0.6,0.8,1.0,1.2,2.0,3.0;
	    
fig = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], xlabel="omn", ylabel="γ/ωᵣ - kyρ = 0.4", linewidth=linewidth)
ax2 = Axis(fig[2,1], xlabel="omn", ylabel="γ/ωᵣ - kyρ = 0.6", linewidth=linewidth)
ax3 = Axis(fig[3,1], xlabel="omn", ylabel="γ/ωᵣ - kyρ = 0.8", linewidth=linewidth)
ax4 = Axis(fig[1,2], xlabel="omn", ylabel="γ/ωᵣ - kyρ = 1.0", linewidth=linewidth)
ax5 = Axis(fig[2,2], xlabel="omn", ylabel="γ/ωᵣ - kyρ = 1.2", linewidth=linewidth)  

lines!(ax1, omn, g1, color = :red,linewidth=linewidth)
lines!(ax2, omn, g2, color = :red,linewidth=linewidth)
lines!(ax3, omn, g3, color = :red,linewidth=linewidth)
lines!(ax4, omn, g4, color = :red,linewidth=linewidth)
lines!(ax5, omn, g5, color = :red,linewidth=linewidth)

lines!(ax1, omn, o1, color = :blue, linewidth=linewidth)
lines!(ax2, omn, o2, color = :blue, linewidth=linewidth)
lines!(ax3, omn, o3, color = :blue, linewidth=linewidth)
lines!(ax4, omn, o4, color = :blue, linewidth=linewidth)
lines!(ax5, omn, o5, color = :blue, linewidth=linewidth)

hlines!(ax1, [0], color = :black, linewidth = linewidth, linestyle = :dash) 
vlines!(ax1, [0], color = :black, linewidth = linewidth, linestyle = :dash)
hlines!(ax2, [0], color = :black, linewidth = linewidth, linestyle = :dash) 
vlines!(ax2, [0], color = :black, linewidth = linewidth, linestyle = :dash)
hlines!(ax3, [0], color = :black, linewidth = linewidth, linestyle = :dash) 
vlines!(ax3, [0], color = :black, linewidth = linewidth, linestyle = :dash)
hlines!(ax4, [0], color = :black, linewidth = linewidth, linestyle = :dash) 
vlines!(ax4, [0], color = :black, linewidth = linewidth, linestyle = :dash)
hlines!(ax5, [0], color = :black, linewidth = linewidth, linestyle = :dash) 
vlines!(ax5, [0], color = :black, linewidth = linewidth, linestyle = :dash)

save("tem_omn_ky_scan.png", fig)

fig




