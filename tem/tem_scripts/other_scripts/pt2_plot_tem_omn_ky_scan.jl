using DelimitedFiles, CairoMakie, GLMakie
GLMakie.activate!()

sc = readdlm("omn_ky_scan")

markersize = 80
linewidth = 5
fontsize = 50

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
	    
fig = Figure()

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax6 = Axis(fig[1,1], xlabel="omn", ylabel="γ/ωᵣ - kyρ = 2.0", linewidth=linewidth)  
ax7 = Axis(fig[2,1], xlabel="omn", ylabel="γ/ωᵣ - kyρ = 3.0", linewidth=linewidth)  

lines!(ax6, omn, g6, color = :red,linewidth=linewidth)
lines!(ax7, omn, g7, color = :red,linewidth=linewidth)

lines!(ax6, omn, o6, color = :blue, linewidth=linewidth)
lines!(ax7, omn, o7, color = :blue, linewidth=linewidth)

hlines!(ax6, [0], color = :black, linewidth = linewidth, linestyle = :dash) 
vlines!(ax6, [0], color = :black, linewidth = linewidth, linestyle = :dash)
hlines!(ax7, [0], color = :black, linewidth = linewidth, linestyle = :dash) 
vlines!(ax7, [0], color = :black, linewidth = linewidth, linestyle = :dash)

save("tem_omn_ky_scan.png", fig)

fig




